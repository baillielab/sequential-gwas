include { get_prefix; get_julia_cmd } from '../modules/utils.nf'

process ImputeGenotypes {
    input:
        path topmed_api_token_file
        path genotypes

    script:
        genotypes_prefix = get_prefix(genotypes[0])
        """${get_julia_cmd} impute ${genotypes_prefix} \
            ${topmed_api_token_file} \
            --password ${params.TOPMED_ENCRYPTION_PASSWORD} \
            --max-concurrent-submissions ${params.TOPMED_MAX_PARALLEL_JOBS} \
            --refresh-rate ${params.TOPMED_REFRESH_RATE} \
            --r2 ${params.IMPUTATION_R2_FILTER} \
            --samples-per-file ${params.N_SAMPLES_PER_IMPUTATION_JOBS}
        """
}


process GetChromosomes {
    input:
        path genotypes

    output:
        path "chromosomes.txt"

    script:
        input_prefix = get_prefix(genotypes[0])
        """
        ${get_julia_cmd(task.cpus)} write-chromosomes ${input_prefix}
        """
}

workflow Imputation {
    topmed_api_token = file("assets/topmed-api-token")
    bed_genotypes = Channel.fromPath("${params.GENOTYPES_PREFIX}.{bed,bim,fam}").collect()
    ImputeGenotypes(topmed_api_token, bed_genotypes)
}