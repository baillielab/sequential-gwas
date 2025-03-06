include { get_prefix } from '../modules/utils.nf'

process MakeGWASGroups {
    publishDir "results/groups", mode: 'symlink'

    input:
        path(covariates)
        path(variables_config)

    output:
        path("group*")

    script:
        """
        ${params.JULIA_CMD} make-gwas-groups \
            ${covariates} \
            ${variables_config} \
            --output-prefix group
        """
}

process RegenieMAFMACSNPs {
    publishDir "results/qced_genotypes", mode: 'symlink'
    input:
        tuple path(genotypes_bed), path(genotypes_bim), path(genotypes_fam)

    output:
        path("qc_pass.snplist")

    script:
        genotypes_prefix = get_prefix(genotypes_bed)
        """
        plink2 \
            --bfile ${genotypes_prefix} \
            --maf ${params.REGENIE_MAF} \
            --mac ${params.REGENIE_MAC} \
            --write-snplist --no-id-header \
            --out qc_pass
        """
}

process RegenieStep1 {
    input:
        tuple path(genotypes_bed), path(genotypes_bim), path(genotypes_fam)
        tuple path(variants)
        path phenotypes
        path covariates
        val phenotypes_type

    output:
        path "regenie_step1_${phenotypes_type}*"

    script:
        genotypes_prefix = get_prefix(genotypes_bed)
        """
        regenie \
            --step 1 \
            --bed ${genotypes_prefix} \
            --extract ${variants} \
            --phenoFile ${phenotypes} \
            --covarFile ${covariates} \
            --${phenotypes_type} \
            --bsize ${params.REGENIE_BSIZE} \
            --lowmem \
            --out regenie_step1_${phenotypes_type}
        """
}

workflow GWAS {
    genotypes = Channel.fromPath("${params.GENOTYPES_PREFIX}*{bed,bim,fam}").collect()
    covariates = file(params.COVARIATES, checkIfExists: true)
    variables_config = file(params.VARIABLES_CONFIG, checkIfExists: true)
    MakeGWASGroups(covariates, variables_config)
    RegenieMAFMACSNPs(genotypes)
    // RegenieStep1(genotypes, RegenieMAFMACSNPs.out, params.VARIANTS, params.PHENOTYPES, params.COVARIATES, "phenoFile", "phenoCol")
}