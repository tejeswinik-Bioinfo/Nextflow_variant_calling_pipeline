#!/usr/bin/env nextflow

process gatk_mutect2 {

    publishDir "${params.variant_call}/${sample_id}", mode: "copy"
    conda "bioconda::gatk4=4.6.2.0"
    
    input:
    tuple val (metadata), path (tumor_bam), path (tumor_bai)

    output:
    tuple val (metadata), path ("*somatic*"), emit: vcf
    tuple val (metadata), path ("*somatic*.tbi"), emit: vcf_index
    path("*.stats"), emit: stats
    tuple val (metadata), path ("*f1r2*"), emit: f1r2

    script:
    sample_id = metadata.sampleName

    """
    gatk Mutect2 \
    -R ${params.ref} \
    -I ${tumor_bam} \
    --germline-resource ${params.gNOMAD} \
    --panel-of-normals ${params.PON} \
    --f1r2-tar-gz ${sample_id}_f1r2.tar.gz \
    -O ${sample_id}_somatic.vcf.gz

    """
}

process gatk_mutect2_tumor_normal {

    publishDir "${params.variant_call}/${sample_id}", mode: "copy"
    conda "bioconda::gatk4=4.6.2.0"
    
    input:
    tuple val (metadata), path (tumor_bam), path (normal_bam)

    output:
    tuple val (metadata), path ("*somatic.vcf.gz"), emit: vcf
    tuple val (metadata), path ("*somatic.vcf.gz.tbi"), emit: vcf_index
    path("*.stats"), emit: stats
    path ("*f1r2*"), emit: f1r2

    script:
    tumor_id = metadata.sampleName
    normal_id = normal_bam.basename

    """
    gatk Mutect2 \
    -R ${params.ref} \
    -I ${tumor_bam} \
    -I ${normal_bam} \
    -normal ${normal_id} \
    --germline-resource ${params.gNOMAD} \
    --panel-of-normals ${params.PON} \
    --f1r2-tar-gz ${tumor_id}_f1r2.tar.gz \
    -O ${tumor_id}_somatic.vcf.gz

    """
}

process gatk_getpileupsummaries{
    publishDir "${params.variant_call}/${sample_id}", mode: "copy"
    conda "bioconda::gatk4=4.6.2.0"

    input:
    tuple val (metadata), path (bam_file), path (bai_file)
    output:
    tuple val (metadata), path ("*pileup_summary*")

    script:
    sample_id = metadata.sampleName
    """
    export JAVA_OPTS="-Xmx59G"
    gatk --java-options "-Xmx59G" GetPileupSummaries \
    --verbosity DEBUG \
    -I ${bam_file} \
    -R ${params.ref} \
    -L ${params.common_variants} \
    -V ${params.common_variants} \
    -O ${sample_id}_pileup_summary.table \
    2>&1 | tee -a gatk_debug_live.log
    """
}

process gatk_calculatecontamination{
    publishDir "${params.variant_call}/${sample_id}", mode: "copy"
    conda "bioconda::gatk4=4.6.2.0"

    input:
    tuple val (metadata), path (pileup_table)
    output:
    tuple val (metadata), path ("*contamination*")

    script:
    sample_id = metadata.sampleName
    """
    gatk CalculateContamination \
    -I ${pileup_table[0]} \
    -O ${sample_id}_contamination.table
    """

}

process gatk_orientationbias{
    publishDir "${params.variant_call}/${sample_id}", mode: "copy"
    conda "bioconda::gatk4=4.6.2.0"

    input:
    tuple val (metadata), path (f1r2_file)
    output:
    tuple val (metadata), path ("*orientation_bias*")

    script:
    sample_id = metadata.sampleName
    """
    gatk LearnReadOrientationModel \
    -I ${f1r2_file} \
    -O ${sample_id}_orientation_bias.tar.gz
    """ 
}

process gatk_filtermutectcalls{
    publishDir "${params.variant_call}/${sample_id}", mode: "copy"
    conda "bioconda::gatk4=4.6.2.0"

    input:
    tuple val (metadata), path (raw_vcf), path (orientation_bias), path (contamination_table)

    output:
    tuple val (metadata), path ("*filtered_variants.vcf.gz"), path ("*filtered_variants.vcf.gz.tbi")

    script:
    sample_id = metadata.sampleName
    """
    gatk FilterMutectCalls \
    -V ${raw_vcf[0]} \
    -R ${params.ref} \
    --contamination-table ${contamination_table} \
    --ob-priors ${orientation_bias} \
    -O ${sample_id}_filtered_variants.vcf.gz
    """ 
}

process gatk_funcotator_datasource_downloader{
    storeDir "${ds_parent}", mode: "copy"
    conda "bioconda::gatk4=4.6.2.0"

    output:
    path ("FUNCOTATOR_DATASOURCE"), emit: ds_dir

    script:
    ds_parent = file(params.funcotator_datasource).getParent()

    """
    gatk FuncotatorDataSourceDownloader \
    --somatic \
    --hg38 \
    --extract-after-download \
    -O FUNCOTATOR_DATASOURCE
    """
}


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
    -V ${filtered_vcf[0]} \
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
    publishDir "${params.outdir}/annotation/${sample_id}", mode: 'copy'
    conda "bioconda::snpeff=5.4.0c bioconda::htslib=1.19"

    input:
    tuple val(cohort_metadata), path(filtered_extracted_vcf)
    path(snpeff_db) // This comes from the download process

    output:
    tuple val(cohort_metadata), path("${filtered_extracted_vcf.simpleName}_ann.vcf.gz"), emit: ann_vcf
    tuple val(cohort_metadata), path("${filtered_extracted_vcf.simpleName}_ann.vcf.gz.tbi"), emit: ann_vcf_index

    script:
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
