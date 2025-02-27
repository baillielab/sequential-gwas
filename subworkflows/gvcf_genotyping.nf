include { IndexReferenceGenome } from '../modules/index_ref_genome.nf'
include { CreateSequenceDictionary } from '../modules/gatk_create_seq_dict.nf'
include { GenotypeGVCFs } from '../modules/genotype_gvcfs.nf'
include { MergeGenotypes } from '../modules/merge_plink_files.nf'
include { get_prefix } from '../modules/utils.nf'

workflow GVCFGenotyping {
    take:
        wgs_gvcfs
        shared_variants_gatk
        shared_variants_plink
        reference_genome
        kgp_bim
        
    main:
        reference_genome_index = IndexReferenceGenome(reference_genome)
        reference_genome_dict = CreateSequenceDictionary(reference_genome)
        genotypes = GenotypeGVCFs(
            wgs_gvcfs, 
            shared_variants_gatk,
            shared_variants_plink,
            reference_genome, 
            reference_genome_index, 
            reference_genome_dict,
            kgp_bim
        )
        merge_list = genotypes
            .map { it -> get_prefix(it[0].getName()) }
            .collectFile(name: "wgs_merge_list.txt", newLine: true)
        merged = MergeGenotypes(
            genotypes.collect(), 
            merge_list, 
            "${params.WGS_PUBLISH_DIR}/merged",
            "wgs.merged"
        )
    
    emit:
        genotypes = merged.genotypes
        counts = merged.stats
}