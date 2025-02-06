include { LiftOver; PedToBed } from '../modules/misc.nf'
include { QCRawGenotypes } from '../modules/qc.nf'
include { MergeGenotypes } from '../modules/merge_plink_files.nf'

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