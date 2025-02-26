include { get_prefix } from './utils.nf'

process LDPruning {
    input:
        tuple path(bed_file), path(bim_file), path(fam_file)
        path high_ld_regions

    output:
        tuple path("${output_prefix}.bed"), path("${output_prefix}.bim"), path("${output_prefix}.fam")

    script:
        input_prefix = get_prefix(bed_file)
        output_prefix = "${input_prefix}.ldpruned"
        """
        plink2 \
            --bfile ${input_prefix} \
            --maf ${params.PCA_MAF} \
            --indep-pairwise ${params.IP_VALUES} \
            --exclude range ${high_ld_regions} \
            --make-bed \
            --out ${output_prefix}
        """
}