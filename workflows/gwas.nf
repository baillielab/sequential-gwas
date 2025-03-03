include { get_prefix } from './modules/utils.nf'


def write_chr_list(mergeList, genotypes) {
    def bed_files_with_chr = []
    for (bed_file in genotypes.findAll() { it.toString().endsWith(".bed") }) {
        bed_file_str = bed_file.toString()
        matcher = bed_file_str =~ /chr[0-9]{1,2}/
        bed_files_with_chr.add([matcher[0][3..-1].toInteger(), bed_file_str.take(bed_file_str.lastIndexOf('.'))])
    }
    bed_files_with_chr.sort{ it[0] }
    bed_files_with_chr.each { it -> mergeList << "${it[1]}.bed ${it[1]}.bim ${it[1]}.fam\n" }
}

process MakeMergeList {
    publishDir "results/merged_genotypes", mode: 'symlink'
    input:
        path genotypes
    output:
        path "merge_list.txt"
    exec:
        def mergeList = task.workDir.resolve("merge_list.txt")
        write_chr_list(mergeList, genotypes)
}

process MergeBEDs {
    publishDir "results/merged_genotypes", mode: 'symlink'
    input:
        path genotypes
        path merge_list

    output:
        tuple path("genotypesAllChr.bed"), path("genotypesAllChr.bim"), path("genotypesAllChr.fam")

    script:
        """
        plink2 --pmerge-list merge_list.txt --make-bed --out genotypesAllChr
        """
}


process RegenieQC {
    publishDir "results/qced_genotypes", mode: 'symlink'
    input:
        tuple path(genotypes_bed), path(genotypes_bim), path(genotypes_fam)

    output:
        tuple path("qc_pass.snplist"), path("qc_pass.id")

    script:
        genotypes_prefix = get_prefix(genotypes_bed)
        """
        plink2 \
            --bfile ${genotypes_prefix} \
            --maf ${params.QC_MAF} \
            --mac ${params.QC_MAC} \
            --geno ${params.QC_GENOTYPE_MISSING_RATE} \
            --hwe ${params.QC_HWE_P} ${params.QC_HWE_K}\
            --mind ${params.QC_INDIVIDUAL_MISSING_RATE} \
            --write-snplist --write-samples --no-id-header \
            --out qc_pass
        """
}

process RegenieStep1 {
    input:
        tuple path(genotypes_bed), path(genotypes_bim), path(genotypes_fam)
        tuple path(qc_snplist), path(qc_individuals)
        path phenotypes
        path covariates
        path qc_snplist
        path qc_individuals
        val phenotypes_type

    output:
        path "regenie_step1_${phenotypes_type}*"

    script:
        genotypes_prefix = get_prefix(genotypes_bed)
        """
        regenie \
            --step 1 \
            --bed ${genotypes_prefix} \
            --extract ${qc_snplist} \
            --keep ${qc_individuals} \
            --phenoFile ${phenotypes} \
            --covarFile ${covariates} \
            --${phenotypes_type} \
            --bsize ${params.REGENIE_BSIZE} \
            --lowmem \
            --out regenie_step1_${phenotypes_type}
        """
}