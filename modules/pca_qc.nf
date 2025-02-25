include { get_prefix } from './utils.nf'

process PCAQC {
    publishDir "${params.MERGED_PUBLISH_DIR}/pca_qced", mode: 'symlink'

    input:
        tuple path(bed_file), path(bim_file), path(fam_file)
        path ancestry

    output:
        tuple path("${output_prefix}.bed"), path("${output_prefix}.bim"), path("${output_prefix}.fam"), emit: genotypes
        path("*.png"), emit: plots
        path "variants_to_exclude.txt", emit: excluded_variants

    script:
        input_prefix = get_prefix(bed_file)
        output_prefix = "${input_prefix}.after_pca_qc"
        """
        ${params.JULIA_CMD} pca-qc \
            --input-prefix ${input_prefix} \
            --output-prefix ${output_prefix} \
            --ancestry-file ${ancestry} \
            --npcs ${params.N_PCS} \
            --iqr-factor ${params.IQR_FACTOR}
        """
}