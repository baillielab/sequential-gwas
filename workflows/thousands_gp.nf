include { DownloadOrAccessStoredResource } from '../modules/download_or_access.nf'
include { CleanTounsandsGP } from '../modules/clean_thousands_gp.nf'
include { VCFToBed } from '../modules/vcf_to_bed.nf'

workflow ThousandsGP {
    chrs = Channel.of(1..22)
    // Question: How did they know these VCF files (instead of the 1000GP release) would contain our variants?
    chr_vcfs = chrs
        .map { it -> [
            "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased", 
            "CCDG_14151_B01_GRM_WGS_2020-08-05_chr${it}.filtered.shapeit2-duohmm-phased.vcf.gz"
            ]}
    chr_vcfs_idx = chrs
        .map { it -> [
            "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased", 
            "CCDG_14151_B01_GRM_WGS_2020-08-05_chr${it}.filtered.shapeit2-duohmm-phased.vcf.gz.tbi"
            ]} 
    // panel = Channel.of(["ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502", "integrated_call_samples_v3.20130502.ALL.panel"])
    ped = Channel.of(["ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage", "20130606_g1k_3202_samples_ped_population.txt"])
    all_1000P_files = chr_vcfs.concat(chr_vcfs_idx).concat(ped)

    resources = DownloadOrAccessStoredResource(all_1000P_files, params.THOUSANDSGP_DIR)
    vcf_files = resources
        .filter{ it.toString().contains("vcf") }
        .map{ [it.getName().replace(".tbi", ""), it] }
        .groupTuple(sort: true, size: 2)
        .map{ it[1] }

    clean_vcf_files = CleanTounsandsGP(vcf_files)
    VCFToBed(clean_vcf_files)
}