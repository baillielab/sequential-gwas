process GenotypeGVCFs{
    publishDir "${params.WGS_PUBLISH_DIR}/shared_variants", mode: 'symlink'

    input:
        tuple val(prefix), path(gvcf)
        path join_variants_gatk
        path joint_variants_bim
        path reference_genome
        path reference_genome_index
        path reference_genome_dict

    output:
        tuple path("${prefix}.shared.bed"), path("${prefix}.shared.bim"), path("${prefix}.shared.fam")

    script:
        """
        gatk GenotypeGVCFs \
            -R ${reference_genome} \
            -V ${prefix}.gvcf.gz \
            --intervals ${join_variants_gatk} \
            --interval-padding 10 \
            --force-output-intervals ${join_variants_gatk} \
            -O ${prefix}.shared.vcf.gz
        /opt/miniforge3/bin/mamba run -n plink2_env plink2 \
            --vcf ${prefix}.shared.vcf.gz \
            --output-chr chr26 \
            --set-all-var-ids @:# \
            --make-bed \
            --out ${prefix}.shared
        """
}