include { GVCFGenotyping} from '../subworkflows/gvcf_genotyping.nf'
include { GenotypesQC } from '../subworkflows/genotyping_arrays_qc.nf'
include { KGP } from './kgp.nf'
include { DownloadOrAccessReferenceGenome } from '../modules/download_reference_genome.nf'
include { MergeGenotypesAndQC; MergeGenotypingArraysAndWGS } from '../subworkflows/merge_genotypes.nf'
include { WGSIndividuals } from '../modules/wgs_individuals.nf'
include { Report } from '../subworkflows/report.nf'
include { CombineCovariates } from '../modules/combine_covariates.nf'

workflow CombineDatasets {
    main:
        // Process 1000GP dataset
        kgp = KGP()
        kgp_genotypes = kgp.genotypes.map{ it -> it[0..2] }
        kgp_bim_afreq = kgp.genotypes.map{ it -> [it[1], it[3]] }
        // Reference Genome
        reference_genome = DownloadOrAccessReferenceGenome()
        // GRC37 Genotypes
        r8_array = Channel.of(
            "release-r8", 
            file("${params.R8_GENOTYPES}.ped", 
            checkIfExists: true), file("${params.R8_GENOTYPES}.map", checkIfExists: true)
            )
            .collect()
        before_2024_array = Channel.of(
            "release-2021-2023", 
            file("${params.BEFORE_2024_GENOTYPES}.ped", checkIfExists: true),
            file("${params.BEFORE_2024_GENOTYPES}.map", checkIfExists: true)
            )
            .collect()
        grc37_genotypes = r8_array.concat(before_2024_array)
        // GRC38 Genotypes
        grc38_genotypes = Channel.of(
            "release-2024-now",
            file("${params.SINCE_2024_GENOTYPES}.ped", checkIfExists: true),
            file("${params.SINCE_2024_GENOTYPES}.map", checkIfExists: true)
            )
            .collect()
        // Chain file for GRCh37 genotypes
        chain_file = file(params.GRC37_TO_GRC38_CHAIN_FILE, checkIfExists: true)
        // Variants to be flipped: TODO (unused atm)
        // High LD regions
        high_ld_regions = file(params.HIGH_LD_REGIONS, checkIfExists: true)
        
        if (params.WGS_GVCFS != "") {
            // WGS GVCFs
            wgs_gvcfs = Channel.fromFilePairs("${params.WGS_GVCFS}*{.gvcf.gz,.gvcf.gz.tbi}", checkIfExists: true)
            wgs_sample_ids = WGSIndividuals(wgs_gvcfs.map{it -> it[1]}.collect())
            // Basic QC for genotyping arrays
            qced_genotypes = GenotypesQC(
                grc37_genotypes, 
                grc38_genotypes, 
                chain_file,
                kgp_bim_afreq,
                wgs_sample_ids
            )
            // Genotyping of WGS data based on genotyping arrays variants
            wgs_data = GVCFGenotyping(
                wgs_gvcfs, 
                qced_genotypes.gatk_shared_variants,
                qced_genotypes.plink_shared_variants,
                reference_genome,
            )
            wgs_genotypes = wgs_data.genotypes
            // Merge All Genotypes
            merge_output = MergeGenotypingArraysAndWGS(
                qced_genotypes.genotypes, 
                wgs_genotypes,
                qced_genotypes.plink_shared_variants,
                kgp_genotypes,
                kgp.pedigree,
                high_ld_regions
            )
        }
        else {
            wgs_sample_ids = file("${projectDir}/assets/NO_WGS_SAMPLES.txt")
            wgs_genotypes = file("${projectDir}/assets/NO_WGS_SAMPLES.txt")
            // Basic QC for genotyping arrays
            qced_genotypes = GenotypesQC(
                grc37_genotypes, 
                grc38_genotypes, 
                chain_file,
                kgp_bim_afreq,
                wgs_sample_ids
            )
            // Merge All Genotypes
            merge_output = MergeGenotypesAndQC(
                qced_genotypes.genotypes.map{ it -> it[1..3] }, 
                qced_genotypes.plink_shared_variants,
                kgp_genotypes,
                kgp.pedigree,
                high_ld_regions
            )
        }
        // Combine Covariates
        combined_covariates = CombineCovariates(
            merge_output.ancestries,
            merge_output.pcs,
            wgs_sample_ids,
            qced_genotypes.genotypes.filter { it -> it[0] == "release-r8" }.map{ it -> it[3] }.first(),
            qced_genotypes.genotypes.filter { it -> it[0] == "release-2021-2023" }.map{ it -> it[3] }.first(),
            qced_genotypes.genotypes.filter { it -> it[0] == "release-2024-now" }.map{ it -> it[3] }.first()
        )
        // Report
        Report(
            qced_genotypes.unlifted,
            qced_genotypes.initial_bed_files,
            qced_genotypes.basic_qc_logs,
            qced_genotypes.kgp_qc_files_r8,
            qced_genotypes.kgp_qc_files_2021_2023,
            qced_genotypes.kgp_qc_files_2024_now,
            qced_genotypes.plink_shared_variants,
            wgs_genotypes,
            merge_output.merge_log,
            merge_output.unrelated_individuals,
            merge_output.qc_merge_log,
            merge_output.pca_plots,
            merge_output.high_loadings_variants,
            merge_output.final_merged,
            combined_covariates
        )

}