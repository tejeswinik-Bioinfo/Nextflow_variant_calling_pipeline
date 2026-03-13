#!/usr/bin/env nextflow

process bwa_index{
    
    publishDir "${ref_parent}", mode: "copy"
    conda "bioconda::bwa=0.7.19"

    input:
    path ref

    output:
    path "${ref}*"

    script:

    ref_parent = file(params.ref).getParent()
    """
    bwa index ${ref}
    """
}



process bwa_alignment{

    publishDir "${params.outdir}/aligned_reads", mode: "copy"
    conda "bioconda::bwa=0.7.19"

    input:
    tuple val(metadata), path(r1), path(r2)

    output:
    tuple val(metadata), path("*_aligned_reads*")

    script:

    sample_id = metadata.sampleName

    """
    bwa mem \
    -t ${params.threads} \
    -R "@RG\\tID:${sample_id}\\tPL:ILLUMINA\\tSM:${sample_id}" \
    ${params.ref} \
    $r1 \
    $r2 > ${sample_id}_aligned_reads.sam
    """

}

