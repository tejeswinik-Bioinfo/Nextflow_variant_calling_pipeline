#!/usr/bin/env nextflow

process gatk_Mark_Duplicates{

    publishDir "${params.outdir}/data", mode: "copy"
    conda "bioconda::gatk4=4.6.2.0"

    input:
    tuple val(metadata), path (aligned_sam)

    output:
    tuple val(metadata), path ("*sorted_dedup*")

    script:

    sample_id = metadata.sampleName

    """
    gatk MarkDuplicatesSpark \
    -I $aligned_sam \
    -O ${sample_id}_sorted_dedup_reads.bam
    """

}


process gatk_base_recalibrator{

    publishDir "${params.data}", mode: "copy"
    conda "bioconda::gatk4=4.6.2.0"

    input:
    tuple val(metadata), path (dedup_bam)
    output:
    tuple val(metadata), path ("*")

    script:
    sample_id = metadata.sampleName

    """
    gatk BaseRecalibrator \
    -I ${dedup_bam[0]}\
    -R ${params.ref} \
    --known-sites ${params.known_sites} \
    -O ${sample_id}_recal_data.table
    """

}