include { ExtractSharedVariantsFromKGP } from '../modules/extract_shared_variants_from_kgp.nf'
include { MergeGenotypes } from '../modules/merge_plink_files.nf'
include { EstimateAncestry } from '../modules/estimate_ancestry.nf'
include { MakePopFile } from '../modules/make_pop_file.nf'
include { LDPruning } from '../modules/ld_pruning.nf'
include { get_prefix } from '../modules/utils.nf'

workflow AncestryEstimation {
    take:
        genotypes
        shared_variants
        kgp_genotypes
        kgp_pedigree
        high_ld_regions

    main:
        kgp_shared_genotypes = ExtractSharedVariantsFromKGP(
            kgp_genotypes,
            shared_variants
        )
        all_genotypes = genotypes.concat(kgp_shared_genotypes)
        merge_list = all_genotypes
            .map { it -> get_prefix(it[0].getName()) }
            .collectFile(name: "merge_list.txt", newLine: true)
        merged = MergeGenotypes(
            all_genotypes.collect(), 
            merge_list,
            "${params.ANCESTRY_PUBLISH_DIR}/merged",
            "genotypes_and_kgp.merged"
            )
        ld_pruned_genotypes = LDPruning(merged, high_ld_regions)
        ancestry = EstimateAncestry(ld_pruned_genotypes, kgp_pedigree)

    emit:
        ancestry.ancestry
}