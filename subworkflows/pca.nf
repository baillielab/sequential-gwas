include { LDPruning } from '../modules/ld_pruning.nf'
include { PCAFindHighLoadings } from '../modules/pca_qc.nf'
include { FilterHighLoadingsVariants } from '../modules/filter_high_loadings_variants.nf'

workflow PCAQC {
    take:
        genotypes
        ancestry
        high_ld_regions
    main:
        ld_pruned_genotypes = LDPruning(genotypes, high_ld_regions)
        pca_qc_output = PCAFindHighLoadings(ld_pruned_genotypes, ancestry)
        final_genotypes = FilterHighLoadingsVariants(genotypes, pca_qc_output.high_loadings_variants)

    emit:
        genotypes = final_genotypes
        plots = pca_qc_output.plots
        high_loadings_variants = pca_qc_output.high_loadings_variants
        pcs = pca_qc_output.pcs
}