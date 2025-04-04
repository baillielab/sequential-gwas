include { get_prefix; get_julia_cmd } from '../modules/utils.nf'

process TOPMedImputation {
    cpus = params.TOPMED_MAX_PARALLEL_JOBS
    // label "bigmem"

    input:
        path topmed_api_token_file
        path genotypes

    output:
        path "*.txt", emit: jobs_files

    script:
        genotypes_prefix = get_prefix(genotypes[0])
        """
        ${get_julia_cmd(task.cpus)} impute \
            ${genotypes_prefix} \
            ${topmed_api_token_file} \
            --password ${params.TOPMED_ENCRYPTION_PASSWORD} \
            --max-concurrent-submissions ${params.TOPMED_MAX_PARALLEL_JOBS} \
            --refresh-rate ${params.TOPMED_REFRESH_RATE} \
            --r2 ${params.IMPUTATION_R2_FILTER} \
            --samples-per-file ${params.N_SAMPLES_PER_IMPUTATION_JOBS} \
            --output-dir .
        """
}

process DownloadJobResults {
    label "hyperthreaded"

    input:
        path topmed_api_token_file
        val jobs_file
    
    output:
        path "qcreport/*", emit: qcreports
        path "statisticDir/*", emit: statistics
        path "local/*dose.vcf.gz", emit: imputed_genotypes
        path "local/*empiricalDose.vcf.gz", emit: imputed_empirical_genotypes
        path "local/*info.gz", emit: info

    script:
        """
        ${get_julia_cmd(task.cpus)} download-topmed-job \
            ${jobs_file} \
            ${topmed_api_token_file} \
            --password ${params.TOPMED_ENCRYPTION_PASSWORD} \
            --refresh-rate ${params.TOPMED_REFRESH_RATE} \
            --output-dir .
        """

    stub:
        """
        mkdir ${job_id} imputed_1-5001/local imputed_1-5001/qcreport imputed_1-5001/statisticDir
        touch imputed_${job_id}/qcreport/qcreport.html
        touch imputed_${job_id}/statisticDir/chunks-excluded.txt
        touch imputed_${job_id}/statisticDir/snps-excluded.txt
        touch imputed_${job_id}/statisticDir/typed-only.txt
        touch imputed_${job_id}/local/chr_1.dose.vcf.gz imputed_1-5001/local/chr_1.empiricalDose.vcf.gz imputed_1-5001/local/chr_1.info.gz
        touch imputed_${job_id}/local/chr_2.dose.vcf.gz imputed_1-5001/local/chr_2.empiricalDose.vcf.gz imputed_1-5001/local/chr_2.info.gz
        """
}

process MergeVCFsByChr {
    label "hyperthreaded"

    input:
        tuple val(chr), path(vcf_files)

    output:
        path("${output}")

    script:
        output = "${chr}.vcf.gz"
        sorted_vcf_files_string = vcf_files
            .findAll { x -> x.getName().endsWith("vcf.gz") }
            .sort{ x -> x.getName().tokenize("-")[1].toInteger() }
            .join("\n")
        """
        echo "${sorted_vcf_files_string}" > merge_list.txt

        mamba run -n bcftools_env bcftools merge --threads ${task.cpus} -o ${output} -O z -l merge_list.txt
        """
}

process IndexVCF {
    label "multithreaded"

    input:
        path vcf_file

    output:
        path "${vcf_file}.tbi"

    script:
        """
        mamba run -n bcftools_env bcftools index --tbi --threads ${task.cpus} ${vcf_file}
        """
}

process VCFToPGEN {
    label "multithreaded"

    input:
        path vcf_file

    output:
        tuple path("${output_prefix}.pgen"), path("${output_prefix}.pvar"), path("${output_prefix}.psam")

    script:
        output_prefix = get_prefix(get_prefix(vcf_file))
        """
        plink2 --vcf ${vcf_file} --make-pgen --threads ${task.cpus} --out ${output_prefix}
        """
}

process MergePGENsByChr {
    label "multithreaded"

    input:
        tuple val(chr), path(pgen_files)

    output:
        tuple path("${chr}.pgen"), path("${chr}.psam"), path("${chr}.pgi")

    script:
        output = "${chr}.merged.vcf.gz"
        sorted_pgen_files_string = pgen_files
            .findAll { x -> x.getName().endsWith("pgen") } // Only pgen files
            .collect { x -> get_prefix(x) } // Only keep prefix as expected by plink2
            .sort { x -> x.tokenize("-")[1].toInteger() } // Sort by sample batch to retain order of original fileset in imputed merged file
            .join("\n")
        """
        echo "${sorted_pgen_files_string}" > merge_list.txt

        plink2 --set-all-var-ids @:#\\\$r,\\\$a --new-id-max-allele-len 81 --pmerge-list merge_list.txt --threads ${task.cpus} --make-pgen --out ${chr}
        """
}

process MergeVCFsToBGENByChr {
    input:
        tuple val(chr), path(vcf_files)

    output:
        tuple val(chr), path("${output_prefix}.bgen"), path("${output_prefix}.sample")

    script:
        output_prefix = "${chr}.merged"
        qctool_input_string = vcf_files
            .sort{ x -> x.getName().tokenize("-")[1].toInteger() }
            .join(" -g ")
        """
        qctool -g ${qctool_input_string} -og ${output_prefix}.bgen -os ${output_prefix}.sample
        """
}

workflow Impute {
    topmed_api_token = file(params.TOPMED_TOKEN_FILE)
    bed_genotypes = Channel.fromPath("${params.GENOTYPES_PREFIX}.{bed,bim,fam}").collect()

    if (params.TOPMED_JOBS_LIST == "NO_TOPMED_JOBS") {
        jobs_files = TOPMedImputation(topmed_api_token, bed_genotypes)
    }
    else {
        jobs_files = Channel.fromList(params.TOPMED_JOBS_LIST)
            .collectFile { it -> [ "${it}.txt", it + '\n' ] }
    }

    jobs_results = DownloadJobResults(topmed_api_token, jobs_files)
    all_dose_vcfs = jobs_results.imputed_genotypes.flatten()

    // Make BGEN Output
    // Idea 1: Merge vcf and then convert
    all_dose_vcfs_indices = IndexVCF(all_dose_vcfs)
    dose_vcfs_and_indices_by_chr = all_dose_vcfs.concat(all_dose_vcfs_indices)
        .map{ it -> [it.getName().tokenize(".")[1], it]}
        .groupTuple()
        .view()
    chr_vcfs = MergeVCFsByChr(dose_vcfs_and_indices_by_chr)
    pgen_files = VCFToPGEN(chr_vcfs)

    // Idea 2: Convert all VCF files to BGEN and then merge: not working since plink2 can't do non-concatenating merge...
    // bgen_files = VCFToPGEN(all_dose_vcfs)
    // bgen_files_by_chr = bgen_files
    //     .map{ it -> [it[0].getName().tokenize(".")[1], it]} // Prepend Chromosome to list
    //     .groupTuple() // Group by chromosome
    //     .map { it -> [it[0], it[1].flatten()]} // flatten all files in a single list per chromosome
    //     .view()
    // MergePGENsByChr(bgen_files_by_chr)

}
