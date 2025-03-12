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

process RegenieStep1 {
    label "multithreaded"
    publishDir "${params.PUBLISH_DIR}/gwas/${group}/regenie_step1", mode: 'symlink'

    input:
        tuple val(group), path(bed), path(bim), path(fam), path(phenotypes), path(covariates)

    output:
        tuple val(group), path("${group}.step1_1.loco"), path("${group}.step1_pred.listrelative")

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
        awk '{sub(".*/", "", \$2); print \$1, \$2}' ${group}.step1_pred.list > ${group}.step1_pred.listrelative
        """
}

process RegenieStep2 {
    label "multithreaded"
    publishDir "${params.PUBLISH_DIR}/gwas/${group}/regenie_step2", mode: 'symlink'

    input:
        tuple path(bgen), path(bgi), path(sample)
        tuple val(group), path(covariates), path(phenotypes), path(individuals), path(step1_loco), path(step1_pred)

    output:
        path "*.regenie"

    script:
        outprefix = "${get_prefix(bgen)}.${group}.step2"
        """
        mamba run -n regenie_env regenie \
            --step 2 \
            --bgen ${bgen} \
            --bgi ${bgi} \
            --sample ${sample} \
            --phenoFile ${phenotypes} \
            --covarFile ${covariates} \
            --keep ${individuals} \
            --bt \
            --firth --approx --pThresh 0.01 \
            --pred ${step1_pred} \
            --bsize ${params.REGENIE_BSIZE} \
            --out ${outprefix}
        """
}

process MakeGWASPlots {
    publishDir "${params.PUBLISH_DIR}/gwas/${group}/plots", mode: 'symlink'

    input:
        tuple val(group), path(gwas_results)

    output:
        path "*.png"

    script:
        """
        ${get_julia_cmd(task.cpus)} make-gwas-plots ${gwas_results}
        """
}

workflow GWAS {
    // Inputs
    genotypes = Channel.fromPath("${params.GENOTYPES_PREFIX}.{bed,bim,fam}").collect()
    bgen_genotypes = Channel.fromPath("${params.BGEN_GENOTYPES_PREFIX}.{bgen,bgen.bgi,sample}").collect(sort: true)
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
    // Extract covariates and pcs for each group
    group_covariates = group_files
        .filter { group, file -> file.getName().contains("covariates") }
    covariates_and_pcs = group_covariates
        .join(group_pcs)
    group_covariates_and_pcs = MergeCovariatesPCs(covariates_and_pcs)
    // Extract phenotypes for each group
    group_phenotypes = group_files
        .filter { group, file -> file.getName().contains("phenotype") }
    // Run Regenie Step 1
    group_phenotypes = group_files
        .filter { group, file -> file.getName().contains("phenotype") }
    group_files_step_1 = group_beds
        .join(group_phenotypes)
        .join(group_covariates_and_pcs)
    group_step1_output = RegenieStep1(group_files_step_1)
    // Run Regenie Step 2
    group_inputs_step2 = group_covariates_and_pcs
        .join(group_phenotypes)
        .join(group_samples)
        .join(group_step1_output)
    group_gwas_results = RegenieStep2(bgen_genotypes, group_inputs_step2)
    // Plot
    MakeGWASPlots(group_gwas_results)

}