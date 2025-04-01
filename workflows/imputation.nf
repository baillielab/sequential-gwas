include { get_prefix; get_julia_cmd } from '../modules/utils.nf'

process TOPMedImputation {
    label "multithreaded"

    input:
        path topmed_api_token_file
        path genotypes
        path topmed_jobs_file

    output:
        path "*/qcreport/qcreport.html", emit: qcreports
        path "*/statisticDir/*.txt", emit: statistics
        path "*/local/*dose.vcf.gz", emit: imputed_genotypes
        path "*/local/*empiricalDose.vcf.gz", emit: imputed_empirical_genotypes
        path "*/local/*info.gz", emit: info

    script:
        genotypes_prefix = get_prefix(genotypes[0])
        jobs_list_option = topmed_jobs_file.getName() == "NO_TOPMED_JOBS" ? "" : "--jobs-file ${topmed_jobs_file}"
        """
        ${get_julia_cmd(task.cpus)} impute \
            ${genotypes_prefix} \
            ${topmed_api_token_file} \
            --password ${params.TOPMED_ENCRYPTION_PASSWORD} \
            --max-concurrent-submissions ${params.TOPMED_MAX_PARALLEL_JOBS} \
            --refresh-rate ${params.TOPMED_REFRESH_RATE} \
            --r2 ${params.IMPUTATION_R2_FILTER} \
            --samples-per-file ${params.N_SAMPLES_PER_IMPUTATION_JOBS} ${jobs_list_option}
        """

    stub:
        """
        mkdir imputed_1-5001 imputed_1-5001/local imputed_1-5001/qcreport imputed_1-5001/statisticDir
        touch imputed_1-5001/qcreport/qcreport.html
        touch imputed_1-5001/statisticDir/chunks-excluded.txt
        touch imputed_1-5001/statisticDir/snps-excluded.txt
        touch imputed_1-5001/statisticDir/typed-only.txt
        touch imputed_1-5001/local/chr_1.dose.vcf.gz imputed_1-5001/local/chr_1.empiricalDose.vcf.gz imputed_1-5001/local/chr_1.info.gz
        touch imputed_1-5001/local/chr_2.dose.vcf.gz imputed_1-5001/local/chr_2.empiricalDose.vcf.gz imputed_1-5001/local/chr_2.info.gz


        mkdir imputed_5001-10000 imputed_5001-10000/local imputed_5001-10000/qcreport imputed_5001-10000/statisticDir
        touch imputed_5001-10000/qcreport/qcreport.html
        touch imputed_5001-10000/statisticDir/chunks-excluded.txt
        touch imputed_5001-10000/statisticDir/snps-excluded.txt
        touch imputed_5001-10000/statisticDir/typed-only.txt
        touch imputed_5001-10000/local/chr_1.dose.vcf.gz imputed_5001-10000/local/chr_1.empiricalDose.vcf.gz imputed_5001-10000/local/chr_1.info.gz
        touch imputed_5001-10000/local/chr_2.dose.vcf.gz imputed_5001-10000/local/chr_2.empiricalDose.vcf.gz imputed_5001-10000/local/chr_2.info.gz
        """
}

process MergeImputedGenotypesByChr {
    label "multithreaded"

    input:
        tuple val(chr), val(samples), path(files, name:"?/*")

    output:
        tuple path("${output_prefix}.bgen"), path("${output_prefix}.bgi"), path("${output_prefix}.sample")

    script:
        output_prefix = "${chr}.merged"
        samples_string = samples.join("\n")
        files_string = files.join("\n")
        """
        echo "${samples_string}" > samples.txt
        echo "${files_string}" > merge_list.txt

        bcftools merge --threads ${task.cpus} -o ${chr}.merged.vcf.gz -O z -l merge_list.txt
        plink2 --vcf ${chr}.merged.vcf.gz --export bgen-1.2 --out ${output_prefix}
        """

    stub:
        output_prefix = "${chr}.merged"
        samples_string = samples.join("\n")
        files_string = files.join("\n")
        """
        echo "${samples_string}" > samples.txt
        echo "${files_string}" > merge_list.txt
        touch ${output_prefix}.bgen ${output_prefix}.bgi ${output_prefix}.sample
        """
}

workflow Impute {
    topmed_api_token = file(params.TOPMED_TOKEN_FILE)
    topmed_jobs_file = file(params.TOPMED_JOBS_LIST)
    bed_genotypes = Channel.fromPath("${params.GENOTYPES_PREFIX}.{bed,bim,fam}").collect()
    imputation_outputs = TOPMedImputation(topmed_api_token, bed_genotypes, topmed_jobs_file)
    imputed_genotypes_by_chr = imputation_outputs.imputed_genotypes
        .flatten()
        .map{ it -> [it.getName().tokenize(".")[0], it.toString().tokenize("/")[-3], it]}
        .groupTuple()

    MergeImputedGenotypesByChr(imputed_genotypes_by_chr)
}
