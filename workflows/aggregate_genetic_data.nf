include { GVCFGenotyping} from '../subworkflows/gvcf_genotyping.nf'
include { GenotypesQC } from '../subworkflows/genotyping_arrays_qc.nf'
include { KGP } from './kgp.nf'
include { DownloadOrAccessReferenceGenome } from '../modules/download_reference_genome.nf'
include { MergeGenotypingArraysAndWGS } from '../subworkflows/merge_genotypes.nf'
include { PCA } from '../subworkflows/pca.nf'

workflow AggregateGeneticData {
    // Process 1000GP dataset
    kgp_genotypes = KGP()
    kgp_bim = kgp_genotypes.map{ it -> it[1] }
    // Reference Genome
    reference_genome = DownloadOrAccessReferenceGenome()
    // WGS GVCFs
    wgs_gvcfs = Channel.fromFilePairs("${params.WGS_GVCFS}*{.gvcf.gz,.gvcf.gz.tbi}", checkIfExists: true)
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
    // Basic QC for genotyping arrays
    qced_genotypes = GenotypesQC(
        grc37_genotypes, 
        grc38_genotypes, 
        variants_to_flip, 
        chain_file,
        kgp_bim
    )
    // Genotyping of WGS data based on genotyping arrays variants
    wgs_shared_genotypes = GVCFGenotyping(
        wgs_gvcfs, 
        qced_genotypes.gatk_shared_variants,
        qced_genotypes.plink_shared_variants,
        reference_genome,
        kgp_bim
    )
    // Merge All Genotypes
    merged_genotypes = MergeGenotypingArraysAndWGS(
        qced_genotypes.genotypes, 
        wgs_shared_genotypes
    )
    // PCA
    PCA(merged_genotypes.genotypes)
    
}