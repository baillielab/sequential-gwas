process GenotypeGVCFs{
    publishDir "${params.WGS_PUBLISH_DIR}/shared_variants", mode: 'symlink'

    input:
        tuple val(prefix), path(gvcf)
        path join_variants_gatk
        path reference_genome
        path reference_genome_index
        path reference_genome_dict
        path kgp_bim

    output:
        tuple path("${output_prefix}.bed"), path("${output_prefix}.bim"), path("${output_prefix}.fam")

    script:
        output_prefix = "${prefix}.shared"
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
        
        julia --project=/opt/sequential-gwas/ /opt/sequential-gwas/bin/seq-gwas.jl \
            complete-bim-with-ref \
            ${output_prefix}.bim \
            ${kgp_bim} \
            --output ${output_prefix}.bim
        """
}