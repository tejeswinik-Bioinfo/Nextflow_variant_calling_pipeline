#!/usr/bin/env nextflow

params.input = ""
include{fastp_trim_reads} from "./modules/fastp"
include{bwa_index; bwa_alignment} from "./modules/bwa"
include{gatk_Mark_Duplicates; gatk_base_recalibrator; gatk_applyBQSR; alignment_metrics; insert_size_metrics} from "./modules/gatk_preprocessing"
include{samtools_faidx; gatk_sequenceDictionary} from "./modules/index_files"
include{gatk_mutect2} from "./modules/gatk_variant_call"


workflow{

    sample_ch = Channel.fromPath(params.input)
    .splitCsv(header: true)
    .map{row ->
        def meta = row.subMap('sampleName', 'pairedEnd')
        def r1 = row.Read_1
        def r2 = row.Read_2
        [meta,r1,r2]
        
    }

    trimmed_reads = fastp_trim_reads(sample_ch)
    //trimmed_reads.trimmed.view()
    //trimmed_reads.report.view()

    ref = file(params.ref)
    bwa_index_files = [file("${params.ref}.amb"), file("${params.ref}.ann"), file("${params.ref}.pac"), file("${params.ref}.bwt"), file("${params.ref}.sa")]
    index_exists = bwa_index_files.every{it.exists()}

    //println(bwa_index_files)
    //println(index_exists)

    if (!index_exists) {
        bwa_index(params.ref)
    }


    bwa_alignment(fastp_trim_reads.out.trimmed).view()

    //aligned_reads = bwa_alignment(trimmed_reads.out.trimmed)
    //aligned_reads.view()

    gatk_Mark_Duplicates(bwa_alignment.out).view()

    other_index_files = [file("${params.ref}.fai"), file("${params.ref}.dict")]
    other_index_exists = other_index_files.every{it.exists()}

    if (!other_index_exists) {
        samtools_faidx(params.ref)
        gatk_sequenceDictionary(params.ref)
    }

    gatk_base_recalibrator(gatk_Mark_Duplicates.out).view()
    gatk_applyBQSR(gatk_Mark_Duplicates.out, gatk_base_recalibrator.out).view()
    alignment_metrics(gatk_applyBQSR.out).view()
    insert_size_metrics(gatk_applyBQSR.out).view()

    //Variant caling using mutect2
    gatk_mutect2(gatk_applyBQSR.out).view()
    


}
