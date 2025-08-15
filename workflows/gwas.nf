include { get_prefix; get_julia_cmd } from '../modules/utils.nf'
include { LOCOPCA } from '../subworkflows/pca.nf'
include { MergeCovariatesPCs } from '../modules/merge_covariates_pcs.nf'

process MakeCovariatesAndGroups {
    publishDir "${params.PUBLISH_DIR}/gwas/groups", mode: 'symlink'

    input:
        path(covariates)
        path(inferred_covariates)
        path(variables_config)

    output:
        path("${output_prefix}.covariates.csv"), emit: covariates
        path("${output_prefix}.phenotypes.csv"), emit: phenotypes
        path("${output_prefix}.covariates_list.txt"), emit: covariates_list
        path("*individuals*"), emit: groups_lists

    script:
        output_prefix = "clean"
        inferred_covariates_option = inferred_covariates.getName() == "NO_INFERRED_COVARIATES" ? "" : " --inferred-covariates ${inferred_covariates}"
        """
        ${get_julia_cmd(task.cpus)} make-gwas-groups \
            ${covariates} \
            ${variables_config} \
            --output-prefix ${output_prefix} \
            --min-group-size ${params.MIN_GROUP_SIZE}${inferred_covariates_option}
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
            --threads ${task.cpus} \
            --memory ${task.memory.toMega().toString()} \
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
    publishDir "${params.PUBLISH_DIR}/gwas/${group}/regenie_step_1", mode: 'symlink'

    input:
        path(phenotypes)
        path(covariates)
        val (covariates_list)
        tuple val(group), path(samples), path(bed), path(bim), path(fam)

    output:
        tuple val(group), path("${group}.step1_1.loco"), path("${group}.step1_pred.listrelative")

    script:
        genotypes_prefix = get_prefix(bed)
        """
        conda run -n regenie_env regenie \
            --step 1 \
            --bed ${genotypes_prefix} \
            --keep ${samples} \
            --phenoFile ${phenotypes} \
            --covarFile ${covariates} \
            --covarColList ${covariates_list.join(',')} \
            --cv ${params.REGENIE_CV_NFOLDS} \
            --bt \
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
        tuple val(group), path(individuals), path(step1_loco), path(step1_pred), val(chr), path(imputed_genotypes)
        path(covariates)
        path(phenotypes)
        val(covariates_list)

    output:
        tuple val(group), path("*.regenie")

    script:
        input_prefix = get_prefix(imputed_genotypes[0])
        outprefix = "${chr}.${group}.step2"
        """
        conda run -n regenie_env regenie \
            --step 2 \
            --pgen ${input_prefix} \
            --phenoFile ${phenotypes} \
            --covarFile ${covariates} \
            --covarColList ${covariates_list.join(',')},CHR${chr.replace('chr', '')}_OUT_PC{1:${params.N_PCS}} \
            --keep ${individuals} \
            --bt \
            --firth --approx --pThresh 0.01 \
            --pred ${step1_pred} \
            --bsize ${params.REGENIE_BSIZE} \
            --out ${outprefix} \
            --threads ${task.cpus}
        """
}

process MergeRegenieResults {
    publishDir "${params.PUBLISH_DIR}/gwas/${group}/results", mode: 'copy'

    input:
        tuple val(group), path(group_results)

    output:
        tuple val(group), path("${output}")

    script:
        output = "${group}.csv"
        """
        ${get_julia_cmd(task.cpus)} merge-regenie-chr-results \
            chr --output=${output}
        """
}

process MakeGWASPlots {
    publishDir "${params.PUBLISH_DIR}/gwas/${group}/plots", mode: 'symlink'

    input:
        tuple val(group), path(group_results)

    output:
        path "${group}.manhattan.png"
        path "${group}.qq.png"

    script:
        """
        ${get_julia_cmd(task.cpus)} gwas-plots \
            ${group_results} \
            ${group} \
            --output-prefix ${group}
        """
}

workflow GWAS {
    // Inputs
    genotypes = Channel.fromPath("${params.GENOTYPES_PREFIX}.{bed,bim,fam}", checkIfExists: true).collect(sort: true)
    imputed_genotypes = Channel.fromPath("${params.IMPUTED_GENOTYPES_PREFIX}*")
        .map { it -> [it.getName().tokenize('.')[0], it] }
        .groupTuple(sort: true)
    chromosomes = imputed_genotypes
        .map { it[0] }
    covariates = file(params.COVARIATES, checkIfExists: true)
    inferred_covariates = file(params.INFERRED_COVARIATES, checkIfExists: true)
    variables_config = file(params.VARIABLES_CONFIG, checkIfExists: true)
    high_ld_regions = file(params.HIGH_LD_REGIONS, checkIfExists: true)
    // Define covariates, phenotypes and groups
    MakeCovariatesAndGroups(covariates, inferred_covariates, variables_config)
    covariates = MakeCovariatesAndGroups.out.covariates
    phenotypes = MakeCovariatesAndGroups.out.phenotypes
    covariates_list = MakeCovariatesAndGroups.out.covariates_list
        .splitText()
        .map { it -> it.trim() }
        .collect()
    group_samples = MakeCovariatesAndGroups
        .out
        .groups_lists
        .flatten()
        .map { it -> [it.getName().tokenize('.')[-2], it] }
    // Extract genotypes for each group
    group_beds = BEDGroupsQCed(genotypes, group_samples)
    group_pcs = LOCOPCA(group_beds, high_ld_regions, chromosomes)
    // Merge PCs and covariates
    pcs = group_pcs.map { it[-1] }.collect()
    covariates_and_pcs = MergeCovariatesPCs(covariates, pcs)
    // Extract phenotypes for each group
    group_individuals_and_beds = group_samples.join(group_beds)
    group_step1_output = RegenieStep1(phenotypes, covariates, covariates_list, group_individuals_and_beds)
    // Run Regenie Step 2
    inputs_step2 = group_samples
        .join(group_step1_output)
        .combine(imputed_genotypes)
    outputs_step2 = RegenieStep2(inputs_step2, covariates_and_pcs, phenotypes, covariates_list)
    // Mrge Regenie results
    group_results = MergeRegenieResults(outputs_step2.groupTuple())
    // Plot
    MakeGWASPlots(group_results)
}