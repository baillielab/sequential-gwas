include { LDPruning } from '../modules/ld_pruning.nf'
include { PlinkPCA } from '../modules/pca.nf'
include { PCAQC } from '../modules/pca_qc.nf'
workflow PCA {
    take:
        genotypes
        ancestry
        high_ld_regions
    main:
        ld_pruned_genotypes = LDPruning(genotypes, high_ld_regions)
        pca_qc_output = PCAQC(ld_pruned_genotypes, ancestry)
    emit:
        pca_qc_output.genotypes
}