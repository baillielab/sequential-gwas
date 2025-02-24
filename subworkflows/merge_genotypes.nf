include { get_prefix } from '../modules/utils.nf'
include { MergeGenotypes } from '../modules/merge_plink_files.nf'
include { KingRelationshipInference } from '../modules/king_relationship_inference.nf'
include { QCMergedGenotypes } from '../modules/qc_merged_genotypes.nf'
include { LDPruning } from '../modules/ld_pruning.nf'

workflow MergeGenotypingArraysAndWGS {
    take:
        qced_array_genotypes
        qced_gws_genotypes

    main:
        all_genotypes = qced_array_genotypes
            .map{ it -> it[1..3] }
            .concat(qced_gws_genotypes)
        merge_list = all_genotypes
            .map { it -> get_prefix(it[0].getName()) }
            .collectFile(name: "merge_list.txt", newLine: true)
        merged_genotypes = MergeGenotypes(
            all_genotypes.collect(), 
            merge_list, 
            "${params.MERGED_PUBLISH_DIR}/merged",
            "genotypes.merged"
        )
        unrelated_individuals = KingRelationshipInference(merged_genotypes)
        qced_merged_genotypes = QCMergedGenotypes(merged_genotypes, unrelated_individuals)

    emit:
        genotypes = qced_merged_genotypes.genotypes
        reports = qced_merged_genotypes.reports
}