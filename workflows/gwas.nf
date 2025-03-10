include { get_prefix; get_julia_cmd } from '../modules/utils.nf'
include { PCA } from '../subworkflows/pca.nf'
include { MergeCovariatesPCs } from '../modules/merge_covariates_pcs.nf'

process MakeCovariatesAndGroups {
    publishDir "${params.PUBLISH_DIR}/gwas/groups", mode: 'symlink'

    input:
        path(covariates)
        path(variables_config)

    output:
        path("group*")

    script:
        """
        ${get_julia_cmd(task.cpus)} make-gwas-groups \
            ${covariates} \
            ${variables_config} \
            --output-prefix group \
            --min-group-size ${params.MIN_GROUP_SIZE}
        """
}

process BEDGroupsQCed {
    label "multithreaded"
    publishDir "${params.PUBLISH_DIR}/gwas/${group_name}/bed", mode: 'symlink'

    input:
        tuple path(genotypes_bed), path(genotypes_bim), path(genotypes_fam)
        tuple val(group_name), path(sample_list)

    output:
        tuple val(group_name), path("${group_name}.bed"), path("${group_name}.bim"), path("${group_name}.fam")

    script:
        genotypes_prefix = get_prefix(genotypes_bed)
        """
        plink2 \
            --bfile ${genotypes_prefix} \
            --keep ${sample_list} \
            --maf ${params.REGENIE_MAF} \
            --mac ${params.REGENIE_MAC} \
            --make-bed \
            --out ${group_name}
        """
}

process MakeBGENGroupAndQC {
    publishDir "${params.PUBLISH_DIR}/gwas/${group_name}/bgen", mode: 'symlink'

    input:
        tuple path(genotypes_bgen), path(genotypes_bgi), path(genotypes_sample)
        tuple val(group_name), path(sample_list)

    output:
        tuple val(group_name), path("${group_name}.bgen"), path("${group_name}.bgi"), path("${group_name}.sample")

    script:
        """
        qctool \
            -g ${genotypes_bgen} \
            -s i${genotypes_sample} \
            -incl-samples ${sample_list} \
            -maf ${params.REGENIE_MAF} \
            -exclude-mac ${params.REGENIE_MAC} \
            -og ${group_name}.bgen \
            -os ${group_name}.sample
        """
}

process RegenieStep1 {
    publishDir "${params.PUBLISH_DIR}/gwas/${group}/regenie_step1", mode: 'symlink'

    input:
        tuple val(group), path(bed), path(bim), path(fam), path(phenotypes), path(covariates)

    output:
        path "${group}.step1*"

    script:
        phenotypes_type = "bt"
        genotypes_prefix = get_prefix(bed)
        """
        mamba run -n regenie_env regenie \
            --step 1 \
            --bed ${genotypes_prefix} \
            --phenoFile ${phenotypes} \
            --covarFile ${covariates} \
            --${phenotypes_type} \
            --bsize ${params.REGENIE_BSIZE} \
            --lowmem \
            --out ${group}.step1
        """
}

workflow GWAS {
    // Inputs
    genotypes = Channel.fromPath("${params.GENOTYPES_PREFIX}.{bed,bim,fam}").collect()
    bgen_genotypes = Channel.fromPath("${params.BGEN_GENOTYPES_PREFIX}.{bgen,bgen.bgi,sample}").collect()
    covariates = file(params.COVARIATES, checkIfExists: true)
    variables_config = file(params.VARIABLES_CONFIG, checkIfExists: true)
    high_ld_regions = file(params.HIGH_LD_REGIONS, checkIfExists: true)
    // Define covariates, phenotypes and groups
    MakeCovariatesAndGroups(covariates, variables_config)
    group_files = MakeCovariatesAndGroups
        .out
        .flatten()
        .map { it -> [it.getName().tokenize('.')[-2], it] }
    group_samples = group_files
        .filter { group, file -> file.getName().contains("individuals") }
    // Extract genotypes for each group
    group_beds = BEDGroupsQCed(genotypes, group_samples)
    group_pcs = PCA(group_beds, high_ld_regions)
    group_covariates = group_files
        .filter { group, file -> file.getName().contains("covariates") }
    covariates_and_pcs = group_covariates
        .join(group_pcs)
    group_covariates_and_pcs = MergeCovariatesPCs(covariates_and_pcs)
    // Run Regenie Step 1
    group_phenotypes = group_files
        .filter { group, file -> file.getName().contains("phenotype") }
    group_files_step_1 = group_beds
        .join(group_phenotypes)
        .join(group_covariates_and_pcs)
    RegenieStep1(group_files_step_1)

    // Association testing
    // MakeBGENGroupAndQC(bgen_genotypes, group_samples)

}