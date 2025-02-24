include { LDPruning } from '../modules/ld_pruning.nf'
include { PlinkPCA } from '../modules/pca.nf'
include { PCAPlots } from '../modules/pca_plots.nf'
workflow PCA {
    take:
        genotypes
        ancestry
        high_ld_regions

    main:
        ld_pruned_genotypes = LDPruning(genotypes, high_ld_regions)
        eigens = PlinkPCA(ld_pruned_genotypes)
        PCAPlots(eigens, ancestry)
}