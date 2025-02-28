include { get_prefix } from './utils.nf'

process FilterHighLoadingsVariants {
    publishDir "${params.PUBLISH_DIR}", mode: 'copy'

    input:
        tuple path(bed_file), path(bim_file), path(fam_file)
        path(high_loadings_variants)

    output:
        tuple path("${output_prefix}.bed"), path("${output_prefix}.bim"), path("${output_prefix}.fam")

    script:
        input_prefix = get_prefix(bed_file)
        output_prefix = "genotypes.arrays_wgs.aggregated"
        """
        plink2 \
            --bfile ${input_prefix} \
            --output-chr chr26 \
            --exclude ${high_loadings_variants} \
            --make-bed \
            --out ${output_prefix}
        """
}