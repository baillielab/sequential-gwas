process MergeGenotypes {
    publishDir "${publish_dir}", mode: 'symlink'

    input:
        path genotype_files
        path merge_list
        val publish_dir
        val output_prefix

    output:
        tuple path("${output_prefix}.bed"), path("${output_prefix}.bim"), path("${output_prefix}.fam")

    script:
        """
        /opt/miniforge3/bin/mamba run -n plink_env plink \
            --output-chr chr26 \
            --merge-list ${merge_list} \
            --make-bed \
            --out ${output_prefix}
        """
}