include { get_prefix } from '../modules/utils.nf'
include { QCRawGenotypes; QCMergedGenotypes } from '../modules/qc.nf'
include { MakeSharedVariantsList; ExtractSharedVariants} from '../modules/extract_shared_variants.nf'
include { LiftOver; PedToBed; MergeGenotypes } from '../modules/misc.nf'

workflow GenotypesQC {
    take: 
        grc37_genotypes
        grc38_genotypes
        variants_to_flip
        chain_file

    main:
        lifted_genotypes = LiftOver(grc37_genotypes, chain_file)
        genotypes_bed = PedToBed(lifted_genotypes.genotypes.concat(grc38_genotypes))
        qced_genotypes = QCRawGenotypes(genotypes_bed, variants_to_flip)
    emit:
        genotypes = qced_genotypes.genotypes
        reports = qced_genotypes.reports.collect()

}

workflow CombineGenotypingArrays {
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
    qced_shared_genotypes = ExtractSharedVariants(qced_genotypes.genotypes, shared_variants)
    // Merge Genotypes
    merge_list = qced_shared_genotypes
        .map { it -> get_prefix(it[0].getName()) }
        .collectFile(name: "merge_list.txt", newLine: true)
    merged_genotypes = MergeGenotypes(qced_shared_genotypes.collect(), merge_list)
    qced_merged_genotypes = QCMergedGenotypes(merged_genotypes)
    // Aggregate Reports
}