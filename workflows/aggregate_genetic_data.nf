include { get_prefix } from '../modules/utils.nf'
include { QCMergedGenotypes } from '../modules/qc.nf'
include { MakeSharedVariantsList; ExtractSharedVariantsFromPLINK } from '../modules/extract_shared_variants.nf'
include { MergeGenotypes} from '../modules/merge_plink_files.nf'
include { GVCFGenotyping} from '../subworkflows/gvcf_genotyping.nf'
include { GenotypesQC } from '../subworkflows/genotyping_arrays_qc.nf'
include { KGP } from './kgp.nf'
include { QCFromKGP } from '../modules/qc_from_kgp.nf'
include { FlipAndExtract } from '../modules/flip_and_extract.nf'

workflow AggregateGeneticData {
    // Process 1000GP dataset
    kgp_genotypes = KGP()
    kgp_bim = kgp_genotypes.map{ it -> it[1] }.first()
    // Reference Genome
    reference_genome = file(params.REFERENCE_GENOME, checkIfExists: true)
    // WGS GVCFs
    wgs_gvcfs = Channel.fromFilePairs("${params.WGS_GVCFS}*{.gvcf.gz,.gvcf.gz.tbi}", checkIfExists: true)
    // GRC37 Genotypes
    r8_array = Channel.of(
        "release-r8", 
        file("${params.R8_GENOTYPES}.ped", 
        checkIfExists: true), file("${params.R8_GENOTYPES}.map", checkIfExists: true)
        )
        .collect()

    before_2024_array = Channel.of(
        "release-2021-2023", 
        file("${params.BEFORE_2024_GENOTYPES}.ped", checkIfExists: true),
        file("${params.BEFORE_2024_GENOTYPES}.map", checkIfExists: true)
        )
        .collect()
    grc37_genotypes = r8_array.concat(before_2024_array)
    // GRC38 Genotypes
    grc38_genotypes = Channel.of(
        "release-2024-now",
        file("${params.SINCE_2024_GENOTYPES}.ped", checkIfExists: true),
        file("${params.SINCE_2024_GENOTYPES}.map", checkIfExists: true)
        )
        .collect()

    // LiftOver and QC Genotypes
    variants_to_flip = file(params.VARIANTS_TO_FLIP_GRC38, checkIfExists: true)
    chain_file = file(params.GRC37_TO_GRC38_CHAIN_FILE, checkIfExists: true)
    qced_genotypes = GenotypesQC(grc37_genotypes, grc38_genotypes, variants_to_flip, chain_file)
    qced_bim_files = qced_genotypes.genotypes.branch{ it ->
        release_r8: it[0] == "release-r8"
                        return [it[0], it[2]]
        release_20212023: it[0] == "release-2021-2023"
                        return [it[0], it[2]]
        release_2024_now: it[0] == "release-2024-now"
                        return [it[0], it[2]]
    }
    kgp_qc_all_files = QCFromKGP(
        qced_bim_files.release_r8,
        qced_bim_files.release_20212023,
        qced_bim_files.release_2024_now,
        kgp_bim
    )
    shared_variants_plink = kgp_qc_all_files.shared_variants_plink.first()
    kgp_qc_release_files = kgp_qc_all_files.release_r8
        .concat(kgp_qc_all_files.release_2021_2023)
        .concat(kgp_qc_all_files.release_2024_now)
    qced_arrays = FlipAndExtract(
        qced_genotypes.genotypes.join(kgp_qc_release_files)
            .map{ release_id, bed, bim, fam, flip, new_bim -> [release_id, bed, new_bim, fam, flip]},
        shared_variants_plink
    )
    
    shared_variants_gatk = kgp_qc_all_files.shared_variants_gatk.first()
    wgs_shared_genotypes = GVCFGenotyping(
        wgs_gvcfs, 
        shared_variants_gatk,
        reference_genome,
        kgp_bim
    )

    // Merge Genotypes
    all_genotypes = qced_arrays
        .map{ it -> it[1..3] }
        .concat(wgs_shared_genotypes)
    merge_list = all_genotypes
        .map { it -> get_prefix(it[0].getName()) }
        .collectFile(name: "merge_list.txt", newLine: true)
    merged_genotypes = MergeGenotypes(all_genotypes.collect(), merge_list, params.ARRAY_GENOTYPES_PUBLISH_DIR)
    QCMergedGenotypes(merged_genotypes)
    
    // Aggregate Reports
}