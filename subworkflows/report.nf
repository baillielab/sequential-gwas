include { MakeReport } from '../modules/report.nf'

workflow Report {
    take:
        unlifted_variants_all_releases
        initial_bed_files_all_releases
        basic_qc_logs_all_releases
        kgp_qc_files_r8
        kgp_qc_files_2021_2023
        kgp_qc_files_2024_now
        shared_variants
        wgs
        merge_log
        unrelated_individuals
        qc_merge_log
        pca_plots
        high_loading_variants
        final_merged_genotypes
        covariates

    main:
        unlifted_variants = unlifted_variants_all_releases.branch{ it ->
            release_r8: it[0] == "release-r8"
                        return it[1]
            release_20212023: it[0] == "release-2021-2023"
                        return it[1]
        }
        initial_bed_files = initial_bed_files_all_releases.branch{ it ->
            release_r8: it[0] == "release-r8"
            release_2021_2023: it[0] == "release-2021-2023"
            release_2024_now: it[0] == "release-2024-now"
        }
        basic_qc_logs = basic_qc_logs_all_releases.branch{ it ->
            release_r8: it[0] == "release-r8"
            release_20212023: it[0] == "release-2021-2023"
            release_2024_now: it[0] == "release-2024-now"
        }
        MakeReport(
            unlifted_variants.release_r8,
            unlifted_variants.release_20212023,
            initial_bed_files.release_r8,
            initial_bed_files.release_2021_2023,
            initial_bed_files.release_2024_now,
            basic_qc_logs.release_r8,
            basic_qc_logs.release_20212023,
            basic_qc_logs.release_2024_now,
            kgp_qc_files_r8,
            kgp_qc_files_2021_2023,
            kgp_qc_files_2024_now,
            shared_variants,
            wgs,
            merge_log,
            unrelated_individuals,
            qc_merge_log,
            pca_plots,
            high_loading_variants,
            final_merged_genotypes,
            covariates
        )

}