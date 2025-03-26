include { get_prefix; get_julia_cmd } from './utils.nf'

process MakeReport {
    publishDir "${params.PUBLISH_DIR}", mode: 'copy'

    input:
        path(unlifted_r8)
        path(unlifted_2021_2023)

        tuple val(release_8_id), path(initial_r8_bed), path(initial_r8_bim), path(initial_r8_fam)
        tuple val(release_2021_2023_id), path(initial_2021_2023_bed), path(initial_2021_2023_bim), path(initial_2021_2023_fam)
        tuple val(release_2024_now_id), path(initial_2024_now_bed), path(initial_2024_now_bim), path(initial_2024_now_fam)

        tuple val(release_r8_report), path(release_r8_qc_logs)
        tuple val(release_2021_2023_report), path(release_2021_2023_qc_logs)
        tuple val(release_2024_now_report), path(release_2024_now_qc_logs)

        tuple val(release_r8_kgp_id), path(release_r8_kgp_samples_to_drop), path(release_r8_kgp_flip), path(release_r8_kgp_new_bim)
        tuple val(release_2021_2023_kgp_id), path(release_2021_2023_kgp_samples_to_drop), path(release_2021_2023_kgp_flip), path(release_2021_2023_kgp_new_bim)
        tuple val(release_2024_now_kgp_id), path(release_2024_now_kgp_samples_to_drop), path(release_2024_now_kgp_flip), path(release_2024_now_kgp_new_bim)

        path(shared_variants)
        path(wgs)

        path(merged_genotypes)
        path(unrelated_individuals)
        path(merged_qced_genotypes)
        path(pca_plots)
        path(high_loading_variants)
        path(final_merged_genotypes)
        path(covariates)
    
    output:
        path("report.md")
        path("*.png")

    script:
        initial_bed_r8_prefix = get_prefix(initial_r8_bed)
        initial_bed_2021_2023_prefix = get_prefix(initial_2021_2023_bed)
        initial_bed_2024_now_prefix = get_prefix(initial_2024_now_bed)
        wgs_prefix = get_prefix(wgs[0])
        merged_genotypes_prefix = get_prefix(merged_genotypes[0])
        merged_qced_genotypes_prefix = get_prefix(merged_qced_genotypes[0])
        pca_plots_prefix = get_prefix(get_prefix(get_prefix(pca_plots[0])))
        final_merged_genotypes_prefix = get_prefix(final_merged_genotypes[0])
        """
        ${get_julia_cmd(task.cpus)} make-report \
            ${unlifted_r8} \
            ${unlifted_2021_2023} \
            ${initial_bed_r8_prefix} \
            ${initial_bed_2021_2023_prefix} \
            ${initial_bed_2024_now_prefix} \
            ${release_r8_qc_logs} \
            ${release_2021_2023_qc_logs} \
            ${release_2024_now_qc_logs} \
            ${release_r8_kgp_samples_to_drop} \
            ${release_2021_2023_kgp_samples_to_drop} \
            ${release_2024_now_kgp_samples_to_drop} \
            ${shared_variants} \
            ${wgs_prefix} \
            ${merged_genotypes_prefix} \
            ${unrelated_individuals} \
            ${merged_qced_genotypes_prefix} \
            ${pca_plots_prefix} \
            ${high_loading_variants} \
            ${final_merged_genotypes_prefix} \
            ${covariates}
        """
}