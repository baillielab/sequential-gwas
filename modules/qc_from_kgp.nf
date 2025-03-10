include { get_prefix; get_julia_cmd } from './utils.nf'

process QCFilesFromKGP {
    publishDir "${params.ARRAY_GENOTYPES_PUBLISH_DIR}/qc_files_from_kgp", mode: 'symlink'

    input:
        tuple val(release_r8_id), path(release_r8_bed), path(release_r8_bim), path(release_r8_fam), path(release_r8_afreq)
        tuple val(release_2021_2023_id), path(release_2021_2023_bed), path(release_2021_2023_bim), path(release_2021_2023_fam), path(release_2021_2023_afreq)
        tuple val(release_2024_now_id), path(release_2024_now_bed), path(release_2024_now_bim), path(release_2024_now_fam), path(release_2024_now_afreq)
        tuple path(kgp_bim), path(kgp_afreq)
        path wgs_sample_ids

    output:
        tuple val(release_r8_id), path("${release_r8_prefix}.samples_to_drop.txt"), path("${release_r8_prefix}.flip.txt"), path("${release_r8_prefix}.new.bim"), emit: release_r8
        tuple val(release_2021_2023_id), path("${release_2021_2023_prefix}.samples_to_drop.txt"), path("${release_2021_2023_prefix}.flip.txt"), path("${release_2021_2023_prefix}.new.bim"), emit: release_2021_2023
        tuple val(release_2024_now_id), path("${release_2024_now_prefix}.samples_to_drop.txt"), path("${release_2024_now_prefix}.flip.txt"), path("${release_2024_now_prefix}.new.bim"), emit: release_2024_now
        path "variants_intersection.txt", emit: shared_variants_plink
        path "variants_intersection.bed", emit: shared_variants_gatk
        path "*summary.csv", emit: summaries

    script:
        release_r8_prefix = get_prefix(release_r8_bed)
        release_2021_2023_prefix = get_prefix(release_2021_2023_bed)
        release_2024_now_prefix = get_prefix(release_2024_now_bed)
        kgp_prefix = get_prefix(kgp_bim)
        """
        
        ${get_julia_cmd(task.cpus)} qc-from-kgp \
            --release-r8 ${release_r8_prefix} \
            --release-2021-2023 ${release_2021_2023_prefix} \
            --release-2024-now ${release_2024_now_prefix} \
            --kgp ${kgp_prefix} \
            --wgs-samples-file ${wgs_sample_ids} \
            --outdir .
        """
}