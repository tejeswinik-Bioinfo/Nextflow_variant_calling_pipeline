#!/usr/bin/env nextflow

process bwa_index{
    
    storeDir "${params.ref_parent}"
    conda "bioconda::bwa=0.7.19"

    input:
    path ref

    output:
    path "${ref}*"

    script:

    """
    bwa index ${ref}
    """
}



process bwa_alignment{

    publishDir "${params.outdir}/aligned_reads/${sample_id}", mode: "copy"
    conda "bioconda::bwa=0.7.19"

    input:
    tuple val(metadata), path(r1), path(r2)

    output:
    tuple val(metadata), path("*_aligned_reads*")

    script:

    sample_id = metadata.sampleName
    // Include sample type (tumor/normal) in read group ID if available
    rg_id = metadata.type ? "${sample_id}_${metadata.type}" : "${sample_id}"

    """
    bwa mem \
    -t ${params.threads} \
    -R "@RG\\tID:${rg_id}\\tPL:ILLUMINA\\tSM:${sample_id}" \
    ${params.ref} \
    $r1 \
    $r2 |  samtools view -bS - > ${sample_id}_aligned_reads.bam
    """

}

