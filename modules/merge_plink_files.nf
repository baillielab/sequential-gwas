process MergeGenotypes {
    label "bigmem"
    label "multithreaded"
    publishDir "${publish_dir}", mode: 'symlink'

    input:
        path genotype_files
        path merge_list
        val publish_dir
        val output_prefix

    output:
        tuple path("${output_prefix}.bed"), path("${output_prefix}.bim"), path("${output_prefix}.fam"), emit: genotypes
        path("${output_prefix}.frq.counts"), emit: stats
        path("${output_prefix}.log"), emit: log

    script:
        """
        plink \
            --output-chr chr26 \
            --merge-list ${merge_list} \
            --make-bed \
            --freq counts \
            --out ${output_prefix}
        """
}