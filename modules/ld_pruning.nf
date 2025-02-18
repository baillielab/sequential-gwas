include { get_prefix } from './utils.nf'

process LDPruning {
    publishDir "${params.PCA_PUBLISH_DIR}/ld_pruned", mode: 'symlink'

    input:
        tuple path(bed_file), path(bim_file), path(fam_file)

    output:
        tuple path("${output_prefix}.bed"), path("${output_prefix}.bim"), path("${output_prefix}.fam")

    script:
        input_prefix = get_prefix(bed_file)
        output_prefix = "${input_prefix}.ldpruned"
        """
        /opt/miniforge3/bin/mamba run -n plink2_env plink2 \
            --bfile ${input_prefix} \
            --maf ${params.PCA_MAF} \
            --indep-pairwise ${params.IP_VALUES} \
            --output-chr chr26 \
            --make-bed \
            --out ${output_prefix}
        """
}