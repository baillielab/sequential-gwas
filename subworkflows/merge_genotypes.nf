include { get_prefix } from '../modules/utils.nf'
include { MergeGenotypes } from '../modules/merge_plink_files.nf'
include { KingRelationshipInference } from '../modules/king_relationship_inference.nf'
include { QCMergedGenotypes } from '../modules/qc_merged_genotypes.nf'
include { LDPruning } from '../modules/ld_pruning.nf'
include { PCA } from '../subworkflows/pca.nf'
include { AncestryEstimation } from '../subworkflows/ancestry.nf'

workflow MergeGenotypingArraysAndWGS {
    take:
        qced_array_genotypes
        qced_gws_genotypes
        plink_shared_variants
        kgp_genotypes
        kgp_pedigree
        high_ld_regions

    main:
        all_genotypes = qced_array_genotypes
            .map{ it -> it[1..3] }
            .concat(qced_gws_genotypes)
        merge_list = all_genotypes
            .map { it -> get_prefix(it[0].getName()) }
            .collectFile(name: "merge_list.txt", newLine: true)
        // Merge
        merged_genotypes_output = MergeGenotypes(
            all_genotypes.collect(), 
            merge_list, 
            "${params.MERGED_PUBLISH_DIR}/merged",
            "genotypes.merged"
        )
        // QC Merged
        unrelated_individuals = KingRelationshipInference(merged_genotypes_output.genotypes)
        qced_merged_genotypes = QCMergedGenotypes(merged_genotypes_output.genotypes, unrelated_individuals)
        // Estimate Ancestry
        ancestry = AncestryEstimation(
            qced_merged_genotypes.genotypes,
            plink_shared_variants,
            kgp_genotypes, 
            kgp_pedigree,
            high_ld_regions
        )
        // PCA
        pca_output = PCA(qced_merged_genotypes.genotypes, ancestry, high_ld_regions)

    emit:
        merged = merged_genotypes_output.genotypes
        qced_merged = qced_merged_genotypes.genotypes
        final_merged = pca_output.genotypes
        unrelated_individuals = unrelated_individuals
        pca_plots = pca_output.plots
        high_loadings_variants = pca_output.high_loadings_variants
}