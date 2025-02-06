process GenotypeGVCFs{
    publishDir "results/wgs/shared_variants", mode: 'symlink'

    input:
        tuple val(prefix), path(gvcf)
        path shared_variants
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
            -L ${shared_variants} \
            -O ${prefix}.shared.vcf.gz
        /opt/miniforge3/bin/mamba run -n plink2_env plink2 \
            --vcf ${prefix}.shared.vcf.gz \
            --set-all-var-ids @:# \
            --make-bed \
            --out ${prefix}.shared
        """
}