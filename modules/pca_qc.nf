include { get_prefix; get_julia_cmd } from './utils.nf'

process PCAFindHighLoadings {
    label "multithreaded"
    publishDir "${params.MERGED_PUBLISH_DIR}/pca_qced", mode: 'symlink'
    publishDir "${params.PUBLISH_DIR}/pca", mode: 'symlink', pattern: "*.{eigenval,eigenvec}"

    input:
        tuple path(bed_file), path(bim_file), path(fam_file)
        path ancestry

    output:
        tuple path("${output_prefix}.bed"), path("${output_prefix}.bim"), path("${output_prefix}.fam"), emit: genotypes
        path("*.png"), emit: plots
        path "variants_to_exclude.txt", emit: high_loadings_variants
        tuple path("${output_prefix}.eigenval"), path("${output_prefix}.eigenvec"), emit: pcs

    script:
        input_prefix = get_prefix(bed_file)
        output_prefix = "${input_prefix}.after_pca_qc"
        """
        ${get_julia_cmd(task.cpus)} pca-qc \
            --input-prefix ${input_prefix} \
            --output-prefix ${output_prefix} \
            --ancestry-file ${ancestry} \
            --npcs ${params.N_PCS} \
            --iqr-factor ${params.IQR_FACTOR}
        """
}