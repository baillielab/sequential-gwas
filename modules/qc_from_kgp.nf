include { get_prefix } from './utils.nf'

process QCFilesFromKGP {
    publishDir "${params.ARRAY_GENOTYPES_PUBLISH_DIR}/qc_files_from_kgp", mode: 'symlink'

    input:
        tuple val(release_r8_id), path(release_r8_bim)
        tuple val(release_2021_2023_id), path(release_2021_2023_bim)
        tuple val(release_2024_now_id), path(release_2024_now_bim)
        path kgp_bim

    output:
        tuple val(release_r8_id), path("${release_r8_prefix}.flip.txt"), path("${release_r8_prefix}.new.bim"), emit: release_r8
        tuple val(release_2021_2023_id), path("${release_2021_2023_prefix}.flip.txt"), path("${release_2021_2023_prefix}.new.bim"), emit: release_2021_2023
        tuple val(release_2024_now_id), path("${release_2024_now_prefix}.flip.txt"), path("${release_2024_now_prefix}.new.bim"), emit: release_2024_now
        path "variants_intersection.txt", emit: shared_variants_plink
        path "variants_intersection.bed", emit: shared_variants_gatk
        path "*summary.csv", emit: summaries

    script:
        release_r8_prefix = get_prefix(release_r8_bim)
        release_2021_2023_prefix = get_prefix(release_2021_2023_bim)
        release_2024_now_prefix = get_prefix(release_2024_now_bim)
        """
        julia --project=/opt/sequential-gwas/ /opt/sequential-gwas/bin/seq-gwas.jl \
            qc-from-kgp \
            --release-r8-bim ${release_r8_bim} \
            --release-2021-2023-bim ${release_2021_2023_bim} \
            --release-2024-now-bim ${release_2024_now_bim} \
            --kgp-bim ${kgp_bim} \
            --outdir .
        """
}