#!/usr/bin/env nextflow

process samtools_faidx{
    publishDir "${ref_parent}", mode: "copy"
    conda "bioconda::samtools=1.23"

    input:
    path ref
    output:
    path "${ref}*"

    script:
    ref_parent = file(params.ref).getParent()

    """
    samtools faidx ${ref}
    """

}

process gatk_sequenceDictionary{
    publishDir "${ref_parent}", mode: "copy"
    conda "bioconda::gatk4=4.6.2.0"

    input:
    path ref
    output:
    path "*.dict"

    script:
    ref_parent = file(params.ref).getParent()

    """
    gatk CreateSequenceDictionary \
    -R ${ref} 
    """

}