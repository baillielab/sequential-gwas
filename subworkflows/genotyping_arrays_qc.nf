include { PedToBed } from '../modules/ped_to_bed.nf'
include { LiftOver } from '../modules/liftover.nf'
include { GenotypingArrayBasicQC } from '../modules/qc_genotyping_array.nf'
include { QCFilesFromKGP } from '../modules/qc_from_kgp.nf'
include { FlipAndExtract } from '../modules/flip_and_extract.nf'

workflow GenotypesQC {
    take: 
        grc37_genotypes
        grc38_genotypes
        variants_to_flip
        chain_file
        kgp_bim_afreq
        wgs_sample_ids

    main:
        // Lift over GRCh37 arrays
        grc37_lifted_bed = LiftOver(grc37_genotypes, chain_file)
        // Convert all files to PLINK bed format
        grch38_bed = PedToBed(grc38_genotypes)
        // Perform basic QC for all files
        genotypes_bed = grc37_lifted_bed.genotypes.concat(grch38_bed)
        qced_genotypes = GenotypingArrayBasicQC(genotypes_bed, variants_to_flip)
        // Generates QC files for each array using the 1000 Genomes Project
        qced_bim_afreq_files = qced_genotypes.genotypes.branch{ it ->
            release_r8: it[0] == "release-r8"
            release_20212023: it[0] == "release-2021-2023"
            release_2024_now: it[0] == "release-2024-now"
        }
        kgp_qc_files = QCFilesFromKGP(
            qced_bim_afreq_files.release_r8,
            qced_bim_afreq_files.release_20212023,
            qced_bim_afreq_files.release_2024_now,
            kgp_bim_afreq,
            wgs_sample_ids
        )
        shared_variants_plink = kgp_qc_files.shared_variants_plink.first()
        kgp_qc_release_files = kgp_qc_files.release_r8
            .concat(kgp_qc_files.release_2021_2023)
            .concat(kgp_qc_files.release_2024_now)
        
        qced_genotypes_with_new_bim_and_flip = qced_genotypes.genotypes
            .map{ it -> it[0..3]} // drop afreq
            .join(kgp_qc_release_files) // join with samples to drop, flip files and kgp bim 
            .map{ release_id, bed, bim, fam, samples_to_drop, flip, new_bim -> [release_id, bed, new_bim, fam, flip, samples_to_drop]} // replace bim with new_bim
        // Flip and extract shared variants
        qced_flipped_genotypes = FlipAndExtract(
            qced_genotypes_with_new_bim_and_flip,
            shared_variants_plink
        )

    emit:
        genotypes = qced_flipped_genotypes.genotypes
        gatk_shared_variants = kgp_qc_files.shared_variants_gatk.first()
        plink_shared_variants = shared_variants_plink
        unlifted = grc37_lifted_bed.unlifted
        initial_bed_files = genotypes_bed
        basic_qc_logs = qced_genotypes.logs
        kgp_qc_files_r8 = kgp_qc_files.release_r8
        kgp_qc_files_2021_2023 = kgp_qc_files.release_2021_2023
        kgp_qc_files_2024_now = kgp_qc_files.release_2024_now

}