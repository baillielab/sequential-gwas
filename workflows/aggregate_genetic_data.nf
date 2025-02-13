include { get_prefix } from '../modules/utils.nf'
include { QCMergedGenotypes } from '../modules/qc.nf'
include { MakeSharedVariantsList; ExtractSharedVariantsFromPLINK } from '../modules/extract_shared_variants.nf'
include { MergeGenotypes} from '../modules/merge_plink_files.nf'
include { GVCFGenotyping} from '../subworkflows/gvcf_genotyping.nf'
include { GenotypesQC } from '../subworkflows/genotyping_arrays_qc.nf'
include { KGP } from './kgp.nf'

workflow AggregateGeneticData {
    // Process 1000GP dataset
    kgp_genotypes = KGP()
    // Reference Genome
    reference_genome = file(params.REFERENCE_GENOME, checkIfExists: true)
    // WGS GVCFs
    wgs_gvcfs = Channel.fromFilePairs("${params.WGS_GVCFS}*{.gvcf.gz,.gvcf.gz.tbi}", checkIfExists: true)
    // GRC37 Genotypes
    r8_array = Channel.fromPath("${params.R8_GENOTYPES}*", checkIfExists: true).collect()
    before_2024_array = Channel.fromPath("${params.BEFORE_2024_GENOTYPES}*", checkIfExists: true).collect()
    grc37_genotypes = r8_array.concat(before_2024_array)
    // GRC38 Genotypes
    grc38_genotypes = Channel.fromPath("${params.SINCE_2024_GENOTYPES}*", checkIfExists: true).collect()

    // LiftOver and QC Genotypes
    variants_to_flip = file(params.VARIANTS_TO_FLIP_GRC38, checkIfExists: true)
    chain_file = file(params.GRC37_TO_GRC38_CHAIN_FILE, checkIfExists: true)
    qced_genotypes = GenotypesQC(grc37_genotypes, grc38_genotypes, variants_to_flip, chain_file)

    // Identify Shared Variants and extract them
    shared_variants = MakeSharedVariantsList(qced_genotypes.genotypes.collect())
    qced_shared_genotypes = ExtractSharedVariantsFromPLINK(qced_genotypes.genotypes, shared_variants.joint_variants_plink)
    wgs_shared_genotypes = GVCFGenotyping(
        wgs_gvcfs, 
        shared_variants.joint_variants_gatk, 
        shared_variants.joint_variants_bim,
        reference_genome
    )

    // Merge Genotypes
    all_genotypes = qced_shared_genotypes.concat(wgs_shared_genotypes)
    merge_list = all_genotypes
        .map { it -> get_prefix(it[0].getName()) }
        .collectFile(name: "merge_list.txt", newLine: true)
    merged_genotypes = MergeGenotypes(all_genotypes.collect(), merge_list, params.ARRAY_GENOTYPES_PUBLISH_DIR)
    QCMergedGenotypes(merged_genotypes)
    
    // Aggregate Reports
}