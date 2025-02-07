process MergeGenotypes {
    publishDir "results/merged", mode: 'symlink'

    input:
        path genotype_files
        path merge_list

    output:
        tuple path("genotypes.merged.bed"), path("genotypes.merged.bim"), path("genotypes.merged.fam")

    script:
        """
        /opt/miniforge3/bin/mamba run -n plink_env plink \
            --output-chr chr26 \
            --merge-list ${merge_list} \
            --make-bed \
            --out genotypes.merged
        """
}