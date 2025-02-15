include { PedToBed } from '../modules/ped_to_bed.nf'
include { LiftOver } from '../modules/liftover.nf'
include { QCGenotypingArray } from '../modules/qc_genotyping_array.nf'
include { MergeGenotypes } from '../modules/merge_plink_files.nf'
include { QCFilesFromKGP } from '../modules/qc_from_kgp.nf'
include { FlipAndExtract } from '../modules/flip_and_extract.nf'

workflow GenotypesQC {
    take: 
        grc37_genotypes
        grc38_genotypes
        variants_to_flip
        chain_file
        kgp_bim

    main:
        // Lift over GRCh37 arrays
        lifted_genotypes = LiftOver(grc37_genotypes, chain_file)
        // Convert all files to PLINK bed format
        genotypes_bed = PedToBed(lifted_genotypes.genotypes.concat(grc38_genotypes))
        // Perform basic QC for all files
        qced_genotypes = QCGenotypingArray(genotypes_bed, variants_to_flip)
        // Generates QC files for each array using the 1000 Genomes Project
        qced_bim_files = qced_genotypes.genotypes.branch{ it ->
            release_r8: it[0] == "release-r8"
                            return [it[0], it[2]]
            release_20212023: it[0] == "release-2021-2023"
                            return [it[0], it[2]]
            release_2024_now: it[0] == "release-2024-now"
                            return [it[0], it[2]]
        }
        kgp_qc_files = QCFilesFromKGP(
            qced_bim_files.release_r8,
            qced_bim_files.release_20212023,
            qced_bim_files.release_2024_now,
            kgp_bim
        )
        shared_variants_plink = kgp_qc_files.shared_variants_plink.first()
        kgp_qc_release_files = kgp_qc_files.release_r8
            .concat(kgp_qc_files.release_2021_2023)
            .concat(kgp_qc_files.release_2024_now)
        
        qced_genotypes_with_new_bim_and_flip = qced_genotypes.genotypes
            .join(kgp_qc_release_files)
            .map{ release_id, bed, bim, fam, flip, new_bim -> [release_id, bed, new_bim, fam, flip]}
        // Flip and extract shared variants
        qced_flipped_genotypes = FlipAndExtract(
            qced_genotypes_with_new_bim_and_flip,
            shared_variants_plink
        )

    emit:
        genotypes = qced_flipped_genotypes
        gatk_shared_variants = kgp_qc_files.shared_variants_gatk.first()

}