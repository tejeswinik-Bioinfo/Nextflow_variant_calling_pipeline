#!/usr/bin/env nextflow


process gatk_funcotator{
    publishDir "${params.outdir}/annotation/${sample_id}", mode: "copy"
    conda "bioconda::gatk4=4.6.2.0"

    input:
    tuple val (metadata), path (filtered_vcf), path (filtered_vcf_index)
    
    output:
    tuple val (metadata), path ("*annotated_variants*")

    script:
    sample_id = metadata.sampleName
    """
    
    gatk Funcotator \
    -V ${filtered_vcf} \
    -R ${params.ref} \
    --ref-version hg38 \
    --data-sources-path ${params.funcotator_datasource} \
    -O ${sample_id}_annotated_variants.vcf.gz \
    --output-file-format VCF
    """
}

process DOWNLOAD_SNPEFF_DB {
    storeDir "${params.snpeff_db_dir}/snpeff_cache" // This ensures the download is saved permanently
    conda "bioconda::snpeff=5.4.0c"

    output:
    path "${params.snpeff_db}", emit: snpeff_db_path

    script:
    """
    snpEff download -v ${params.snpeff_db} -dataDir \$(pwd)
    """
}

process SNPEFF_ANNOTATE {
    tag "${filtered_extracted_vcf.simpleName}"
    publishDir "${params.outdir}/annotation/all/${sample_id}", mode: 'copy'
    conda "bioconda::snpeff=5.4.0c bioconda::htslib=1.19"

    input:
    tuple val(metadata), path(filtered_extracted_vcf), path(filtered_extracted_vcf_index)

    output:
    tuple val(metadata), path("${filtered_extracted_vcf.simpleName}_ann.vcf.gz"), path("${filtered_extracted_vcf.simpleName}_ann.vcf.gz.tbi"), emit: ann_vcf
    tuple val (metadata), path("${filtered_extracted_vcf.simpleName}_*"), emit: ann_vcf_report

    script:
    sample_id = metadata.sampleName
    """
    # -dataDir . tells snpEff to look in the current working directory 
    # (where Nextflow linked the db_dir)
    snpEff \
    -Xmx8g \
    -dataDir \$(pwd) \
    ${params.snpeff_db} \
    ${filtered_extracted_vcf} | bgzip > ${filtered_extracted_vcf.simpleName}_ann.vcf.gz
    tabix -p vcf ${filtered_extracted_vcf.simpleName}_ann.vcf.gz
    """
}

process gatk_select_variants_SNPs {

    publishDir "${params.outdir}/annotation/SNPs/${sample_id}/", mode: "copy"
    conda "bioconda::gatk4=4.6.2.0"

    input:
    tuple val (metadata), path (ann_vcf), path (ann_vcf_index)
    
    output:
    tuple val (metadata), path ("*_snp_ann.vcf.gz"), path ("*_snp_ann.vcf.gz.tbi")

    script:
    sample_id = metadata.sampleName

    """
    gatk SelectVariants \
    -R ${params.ref} \
    -V ${ann_vcf} \
    --select-type-to-include SNP \
    -O ${sample_id}_snp_ann.vcf.gz
    """

}

process gatk_select_variants_INDELs {

    publishDir "${params.outdir}/annotation/INDELs/${sample_id}/", mode: "copy"
    conda "bioconda::gatk4=4.6.2.0"

    input:
    tuple val (metadata), path (ann_vcf), path (ann_vcf_index)
    output:
    tuple val (metadata), path ("*_indels_ann.vcf.gz"), path ("*_indels_ann.vcf.gz.tbi")
    script:
    sample_id = metadata.sampleName
    """
    gatk SelectVariants \
    -R ${params.ref} \
    -V ${ann_vcf} \
    --select-type-to-include INDEL \
    -O ${sample_id}_indels_ann.vcf.gz
    """

}