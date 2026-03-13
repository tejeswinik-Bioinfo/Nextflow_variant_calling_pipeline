#!/usr/bin/env nextflow

process gatk_mutect2{

    publishDir "${params.variant_call}", mode: "copy"
    conda "bioconda::gatk4=4.6.2.0"
    
    input:
    tuple val (metadata), path (bam_file)

    output:
    tuple val (metadata), path ("*raw_variants*")

    script:
    sample_id = metadata.sampleName

    """
    gatk Mutect2 \
    -R ${params.ref} \
    -I ${bam_file[1]} \
    --germline-resource ${params.gNOMAD} \
    --panel-of-normals ${params.PON} \
    -O ${sample_id}_raw_variants.vcf.gz

    """
}