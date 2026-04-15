#!/usr/bin/env nextflow

params.input = ""
params.tumorsampleName = ""
params.normalSampleName = ""
params.tumor_R1 = ""
params.tumor_R2 = ""
params.normal_R1 = ""
params.normal_R2 = ""

include {fastp_trim_reads} from "./modules/fastp"
include {bwa_index; bwa_alignment} from "./modules/bwa"
include {gatk_Mark_Duplicates; gatk_base_recalibrator; gatk_applyBQSR; alignment_metrics; insert_size_metrics} from "./modules/gatk_preprocessing"
include {samtools_faidx; gatk_sequenceDictionary} from "./modules/index_files"
include {gatk_mutect2; gatk_mutect2_tumor_normal; gatk_getpileupsummaries; gatk_calculatecontamination; gatk_orientationbias; gatk_filtermutectcalls; gatk_funcotator_datasource_downloader; gatk_funcotator} from "./modules/gatk_variant_call"

workflow {
    is_csv = params.input
    is_tumor_only = params.tumor_R1?.trim() && params.tumor_R2?.trim() && !params.normal_R1?.trim()
    is_tumor_normal = params.tumor_R1?.trim() && params.tumor_R2?.trim() && params.normal_R1?.trim() && params.normal_R2?.trim()

    // ========== CSV Input Mode ==========
    if (is_csv) {
        println "CSV input mode: ${params.input}"
        csv_ch = Channel.fromPath(params.input)
        .splitCsv(header: true)
        .map { row ->
            meta = row.subMap(['sampleName', 'pairedEnd'])
            if (row.type) meta.type = row.type
            tuple(meta, file(row.Read_1), file(row.Read_2))
        }
        csv_ch.view()

        // Determine if tumor-normal mode based on presence of 'normal' type
        is_csv_tumor_normal = csv_ch.toList()
            .collect()
            map {list -> list.any { it[0].type == 'normal' }}
       
        preproc_out = preprocessing(csv_ch)
        variant_calling(preproc_out, is_csv_tumor_normal)
    }
    // ========== CLI Tumor-Only Mode ==========
    else if (is_tumor_only) {
        println "CLI tumor-only mode: ${params.tumorsampleName}"
        tumor_meta = [sampleName: params.tumorsampleName ?: "tumor_sample", pairedEnd: true]
        tumor_ch = Channel.of(tuple(tumor_meta, file(params.tumor_R1), file(params.tumor_R2)))

        preproc_out = preprocessing(tumor_ch)
        variant_calling(preproc_out, false)
    }
    // ========== CLI Tumor-Normal Mode ==========
    else if (is_tumor_normal) {
        println "CLI tumor-normal mode: tumor=${params.tumorsampleName}, normal=${params.normalSampleName}"
        tumor_meta = [sampleName: params.tumorsampleName ?: "tumor_sample", pairedEnd: true]
        tumor_ch = Channel.of(tuple(tumor_meta, file(params.tumor_R1), file(params.tumor_R2)))

        normal_meta = [sampleName: params.normalSampleName ?: "normal_sample", pairedEnd: true]
        normal_ch = Channel.of(tuple(normal_meta, file(params.normal_R1), file(params.normal_R2)))

        combined_ch = tumor_ch.concat(normal_ch)
        preproc_out = preprocessing(combined_ch)
        variant_calling(preproc_out, true)
    }
    else {
        error "Please provide either:\n" +
              "  1. CSV: --input <samples.csv> (supports multiple samples, tumor-only or tumor-normal with 'type' column)\n" +
              "  2. Tumor-only: --sampleName <name> --tumor_R1 <R1.fastq.gz> --tumor_R2 <R2.fastq.gz>\n" +
              "  3. Tumor-Normal: --sampleName <name> --normalSampleName <name> --tumor_R1 <R1.fastq.gz> --tumor_R2 <R2.fastq.gz> --normal_R1 <R1.fastq.gz> --normal_R2 <R2.fastq.gz>"
    }
}

workflow preprocessing {
    take:
        sample_ch

    main:

        
        trimmed_reads = fastp_trim_reads(sample_ch)
        bwa_index(params.ref)
        samtools_faidx(params.ref)
        gatk_sequenceDictionary(params.ref)

        aligned = bwa_alignment(fastp_trim_reads.out.trimmed)
        dedup = gatk_Mark_Duplicates(aligned)
        recal = gatk_base_recalibrator(gatk_Mark_Duplicates.out.dedup_bam)
        bqsr = gatk_applyBQSR(gatk_Mark_Duplicates.out.dedup_bam, recal)
        bqsr.view()
        
        //alignment_metrics(gatk_applyBQSR.out.dedup_bqsr)
        //insert_size_metrics(gatk_applyBQSR.out.dedup_bqsr)

    emit:
        bqsr
}

workflow variant_calling {
    take:
        bqsr_ch
        is_tumor_normal

    main:

        vcf_ch = Channel.empty()

        if (is_tumor_normal) {
            // Split tumor and normal channels
            tumor_sample = bqsr_ch.filter { meta, bam, bai -> meta.type == 'tumor' || (!meta.type && meta.sampleName.contains("tumor")) }
            normal_sample = bqsr_ch.filter { meta, bam, bai -> meta.type == 'normal' || (!meta.type && meta.sampleName.contains("normal")) }

            // Combine tumor and normal into single tuple
            combined = tumor_sample.combine(normal_sample)
                .map { tumor_tuple, normal_tuple ->
                    def tumor_meta = tumor_tuple[0]
                    def tumor_bam = tumor_tuple[1]
                    def normal_bam = normal_tuple[1]
                    tuple(tumor_meta, tumor_bam, normal_bam)
                }

            mutect2_out = gatk_mutect2_tumor_normal(combined)
            vcf_ch = vcf_ch.mix(gatk_mutect2_tumor_normal.out.vcf)

            pileup_out = gatk_getpileupsummaries(tumor_sample)
        } else {
            mutect2_out = gatk_mutect2(bqsr_ch)
            vcf_ch = vcf_ch.mix(gatk_mutect2.out.vcf)
            pileup_out = gatk_getpileupsummaries(bqsr_ch)
        }

    

        contamination_out = gatk_calculatecontamination(pileup_out)
        filtered_out = gatk_filtermutectcalls(vcf_ch, contamination_out)
        gatk_funcotator(filtered_out).view()
}




