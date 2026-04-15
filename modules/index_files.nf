#!/usr/bin/env nextflow

process samtools_faidx{
    storeDir "${params.ref_parent}"
    conda "bioconda::samtools=1.23"

    input:
    path ref
    output:
    path "${ref}*"

    script:

    """
    samtools faidx ${ref}
    """

}

process gatk_sequenceDictionary{
    storeDir "${params.ref_parent}"
    conda "bioconda::gatk4=4.6.2.0"

    input:
    path ref
    output:
    path "*.dict"

    script:

    """
    gatk CreateSequenceDictionary \
    -R ${ref} 
    """

}