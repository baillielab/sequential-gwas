process FlipAndExtract {
    publishDir "${params.ARRAY_GENOTYPES_PUBLISH_DIR}/flipped_and_shared", mode: 'symlink'

    input:
        tuple val(release_id), path(bed_file), path(bim_file), path(fam_file), path(to_flip), path(samples_to_drop)
        path shared_variants

    output:
        tuple val(release_id), path("${output_prefix}.bed"), path("${output_prefix}.bim"), path("${output_prefix}.fam")

    script:
        output_prefix = "${release_id}.flipped.shared"
        """
        plink \
            --bed ${bed_file} \
            --bim ${bim_file} \
            --fam ${fam_file} \
            --flip ${to_flip} \
            --output-chr chr26 \
            --remove ${samples_to_drop} \
            --extract ${shared_variants} \
            --make-bed \
            --out ${output_prefix}
        """
}