#!/usr/bin/env nextflow

process trimgalore_trim_reads{

    publishDir "${params.outdir}/trimmed_reads/${sample_id}", mode: "copy"
    conda "bioconda::trim-galore=0.6.10"

    input:
    tuple val(metadata), path(r1), path(r2)

    output:
    tuple val(metadata), path ("*_val_1.fq.gz"), path ("*_val_2.fq.gz"), emit: "trimmed"
    tuple val(metadata), path("*_report.txt"), emit: "report"

    script:
    
    sample_id = metadata.sampleName

    """
    trim_galore \
    --paired \
    --quality 30 \
    --length 25 \
    --cores ${params.threads} \
    --gzip \
    "${r1}" "${r2}"

    """


}