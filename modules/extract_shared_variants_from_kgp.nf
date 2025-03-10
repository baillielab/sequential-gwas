include { get_prefix } from './utils.nf'

process ExtractSharedVariantsFromKGP {
    label "multithreaded"
    publishDir "${params.ANCESTRY_PUBLISH_DIR}/kgp_shared", mode: 'symlink'

    input:
        tuple path(bed_file), path(bim_file), path(fam_file)
        path shared_variants

    output:
        tuple path("${output_prefix}.bed"), path("${output_prefix}.bim"), path("${output_prefix}.fam")

    script:
        input_prefix = get_prefix(bed_file)
        output_prefix = "${input_prefix}.shared"
        """
        plink2 \
            --bfile ${input_prefix} \
            --output-chr chr26 \
            --extract ${shared_variants} \
            --make-bed \
            --out ${output_prefix}
        """
}