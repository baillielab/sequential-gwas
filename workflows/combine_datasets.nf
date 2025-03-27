include { GVCFGenotyping} from '../subworkflows/gvcf_genotyping.nf'
include { GenotypesQC } from '../subworkflows/genotyping_arrays_qc.nf'
include { KGP } from './kgp.nf'
include { DownloadOrAccessReferenceGenome } from '../modules/download_reference_genome.nf'
include { MergeGenotypesAndQC } from '../subworkflows/merge_genotypes.nf'
include { WGSIndividuals } from '../modules/wgs_individuals.nf'
include { Report } from '../subworkflows/report.nf'
include { CombineCovariates } from '../modules/combine_covariates.nf'

workflow CombineGeneticDatasets {
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
            genotyping_qc_outputs = GenotypesQC(
                grc37_genotypes, 
                grc38_genotypes, 
                chain_file,
                kgp_bim_afreq,
                wgs_sample_ids
            )
            // Genotyping of WGS data based on genotyping arrays variants
            wgs_data = GVCFGenotyping(
                wgs_gvcfs, 
                genotyping_qc_outputs.gatk_shared_variants,
                genotyping_qc_outputs.plink_shared_variants,
                reference_genome,
            )
            wgs_genotypes = wgs_data.genotypes
            qced_genotypes = genotyping_qc_outputs.genotypes
                .map{ it -> it[1..3] }
                .concat(wgs_genotypes)
        }
        else {
            wgs_sample_ids = file("${projectDir}/assets/NO_WGS_SAMPLES.txt")
            wgs_genotypes = file("${projectDir}/assets/NO_WGS_SAMPLES.txt")
            // Basic QC for genotyping arrays
            genotyping_qc_outputs = GenotypesQC(
                grc37_genotypes, 
                grc38_genotypes, 
                chain_file,
                kgp_bim_afreq,
                wgs_sample_ids
            )
            qced_genotypes = genotyping_qc_outputs.genotypes.map{ it -> it[1..3] }
        }
        // Merge All Genotypes
        merge_qc_outputs = MergeGenotypesAndQC(
            qced_genotypes, 
            genotyping_qc_outputs.plink_shared_variants,
            kgp_genotypes,
            kgp.pedigree,
            high_ld_regions
        )
        // Combine Covariates
        combined_covariates = CombineCovariates(
            merge_qc_outputs.ancestries,
            merge_qc_outputs.pcs,
            wgs_sample_ids,
            genotyping_qc_outputs.genotypes.filter { it -> it[0] == "release-r8" }.map{ it -> it[3] }.first(),
            genotyping_qc_outputs.genotypes.filter { it -> it[0] == "release-2021-2023" }.map{ it -> it[3] }.first(),
            genotyping_qc_outputs.genotypes.filter { it -> it[0] == "release-2024-now" }.map{ it -> it[3] }.first()
        )
        // Report
        Report(
            genotyping_qc_outputs.unlifted,
            genotyping_qc_outputs.initial_bed_files,
            genotyping_qc_outputs.basic_qc_logs,
            genotyping_qc_outputs.kgp_qc_files_r8,
            genotyping_qc_outputs.kgp_qc_files_2021_2023,
            genotyping_qc_outputs.kgp_qc_files_2024_now,
            genotyping_qc_outputs.plink_shared_variants,
            wgs_genotypes,
            merge_qc_outputs.merge_log,
            merge_qc_outputs.unrelated_individuals,
            merge_qc_outputs.qc_merge_log,
            merge_qc_outputs.pca_plots,
            merge_qc_outputs.high_loadings_variants,
            merge_qc_outputs.final_merged,
            combined_covariates
        )

}