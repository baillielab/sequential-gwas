include { GVCFGenotyping} from '../subworkflows/gvcf_genotyping.nf'
include { GenotypesQC } from '../subworkflows/genotyping_arrays_qc.nf'
include { KGP } from './kgp.nf'
include { DownloadOrAccessReferenceGenome } from '../modules/download_reference_genome.nf'
include { MergeGenotypingArraysAndWGS } from '../subworkflows/merge_genotypes.nf'
include { WGSIndividuals } from '../modules/wgs_individuals.nf'
include { Report } from '../subworkflows/report.nf'
include { CombineCovariates } from '../modules/combine_covariates.nf'

workflow AggregateGeneticData {
    // Process 1000GP dataset
    kgp = KGP()
    kgp_genotypes = kgp.genotypes.map{ it -> it[0..2] }
    kgp_bim_afreq = kgp.genotypes.map{ it -> [it[1], it[3]] }
    kgp_bim = kgp.genotypes.map{ it -> it[1] }
    // Reference Genome
    reference_genome = DownloadOrAccessReferenceGenome()
    // WGS GVCFs
    wgs_gvcfs = Channel.fromFilePairs("${params.WGS_GVCFS}*{.gvcf.gz,.gvcf.gz.tbi}", checkIfExists: true)
    wgs_sample_ids = WGSIndividuals(wgs_gvcfs.map{it -> it[1]}.collect())
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
    variants_to_flip = file(params.VARIANTS_TO_FLIP_GRC38, checkIfExists: true)
    // High LD regions
    high_ld_regions = file(params.HIGH_LD_REGIONS, checkIfExists: true)
    // Basic QC for genotyping arrays
    qced_genotypes = GenotypesQC(
        grc37_genotypes, 
        grc38_genotypes, 
        variants_to_flip, 
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
        kgp_bim
    )
    // Merge All Genotypes
    merge_output = MergeGenotypingArraysAndWGS(
        qced_genotypes.genotypes, 
        wgs_data.genotypes,
        qced_genotypes.plink_shared_variants,
        kgp_genotypes,
        kgp.pedigree,
        high_ld_regions
    )
    // Combine Covariates
    covariates = file(params.COVARIATES, checkIfExists: true)
    CombineCovariates(
        covariates,
        merge_output.ancestries,
        merge_output.pcs
    )
    // Report
    Report(
        qced_genotypes.unlifted,
        qced_genotypes.initial_bed_files,
        qced_genotypes.basic_qc_reports,
        qced_genotypes.kgp_qc_files_r8,
        qced_genotypes.kgp_qc_files_2021_2023,
        qced_genotypes.kgp_qc_files_2024_now,
        qced_genotypes.plink_shared_variants,
        wgs_data.genotypes,
        merge_output.merged,
        merge_output.unrelated_individuals,
        merge_output.qced_merged,
        merge_output.pca_plots,
        merge_output.high_loadings_variants,
        merge_output.final_merged
    )
}