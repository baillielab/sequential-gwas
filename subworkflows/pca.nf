include { LDPruning; GroupLDPruning } from '../modules/ld_pruning.nf'
include { PCAFindHighLoadings } from '../modules/pca_qc.nf'
include { MaybeFilterHighLoadingsVariants } from '../modules/filter_high_loadings_variants.nf'
include { LocoGroupPCA } from '../modules/pca.nf'

workflow PCAAnalysis {
    take:
        genotypes
        ancestry
        high_ld_regions

    main:
        ld_pruned_genotypes = LDPruning(genotypes, high_ld_regions)
        pca_qc_output = PCAFindHighLoadings(ld_pruned_genotypes, ancestry)
        final_genotypes = MaybeFilterHighLoadingsVariants(genotypes, pca_qc_output.high_loadings_variants)

    emit:
        genotypes = final_genotypes
        plots = pca_qc_output.plots
        high_loadings_variants = pca_qc_output.high_loadings_variants
        pcs = pca_qc_output.pcs
}

// workflow PCA {
//     take:
//         genotypes
//         high_ld_regions

//     main:
//         ld_pruned_genotypes = GroupLDPruning(genotypes, high_ld_regions)
//         pcs = GroupPCA(ld_pruned_genotypes)

//     emit:
//         pcs
// }

workflow LOCOPCA {
    take:
        genotypes
        high_ld_regions
        chromosomes

    main:
        ld_pruned_genotypes = GroupLDPruning(genotypes, high_ld_regions)
        chrs_and_genotypes = chromosomes.combine(ld_pruned_genotypes)
        chrs_and_pcs = LocoGroupPCA(chrs_and_genotypes)

    emit:
        chrs_and_pcs
}