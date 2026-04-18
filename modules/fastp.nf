#!/usr/bin/env nextflow

process fastp_trim_reads{

    publishDir "${params.outdir}/trimmed_reads/${sample_id}", mode: "copy"
    conda "bioconda::fastp=1.1.0"

    input:
    tuple val(metadata), path(r1), path(r2)

    output:
    tuple val(metadata), path("*trimmed_R1.fastq.gz"), path("*trimmed_R2.fastq.gz"), emit: "trimmed"
    path("*_report.html"), emit: "report"

    script:
    
    sample_id = metadata.sampleName

    """
    fastp \
    -i $r1 \
    -I $r2 \
    -o ${sample_id}_trimmed_R1.fastq.gz \
    -O ${sample_id}_trimmed_R2.fastq.gz \
    --detect_adapter_for_pe \
    --html ${sample_id}_report.html


    """





}