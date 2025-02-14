process FlipAndExtract {
    publishDir "${params.ARRAY_GENOTYPES_PUBLISH_DIR}/flipped_and_shared", mode: 'symlink'

    input:
        tuple val(release_id), path(bed_file), path(bim_file), path(fam_file), path(to_flip)
        path shared_variants

    output:
        tuple val(release_id), path("${output_prefix}.bed"), path("${output_prefix}.bim"), path("${output_prefix}.fam")

    script:
        output_prefix = "${release_id}.flipped.shared"
        """
        /opt/miniforge3/bin/mamba run -n plink_env plink \
            --bed ${bed_file} \
            --bim ${bim_file} \
            --fam ${fam_file} \
            --flip ${to_flip} \
            --extract ${shared_variants} \
            --make-bed \
            --out ${output_prefix}
        """

}