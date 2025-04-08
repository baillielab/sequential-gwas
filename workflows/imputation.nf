include { get_prefix; get_julia_cmd } from '../modules/utils.nf'

process WriteImputationSplitLists {
    input:
        path genotypes

    output:
        path "genomicc.chromosomes.txt", emit: chromosomes
        path "*.keep", emit: samples

    script:
        genotypes_prefix = get_prefix(genotypes[0])
        """
        ${get_julia_cmd(task.cpus)} write-imputation-split-lists \
            ${genotypes_prefix} \
            --output-prefix genomicc \
            --n-samples-per-file ${params.N_SAMPLES_PER_IMPUTATION_JOBS}
        """
}

process MakeVCFSplit {
    label "multithreaded"

    input:
        path genotypes
        tuple val(chr), path(samples)

    output:
        path("${output_prefix}.vcf.gz")

    script:
        genotypes_prefix = get_prefix(genotypes[0])
        samples_id = samples.getName().tokenize(".")[1]
        output_prefix = "genomicc.${chr}.${samples_id}"
        """
        plink2 \
            --bfile ${genotypes_prefix} \
            --keep ${samples} \
            --chr ${chr} \
            --export vcf-4.2 id-delim=@ \
            --out ${output_prefix} \
            --threads ${task.cpus} \
            --output-chr chr26
        
        bgzip ${output_prefix}.vcf
        """
}

process TOPMedImputation {
    cpus = params.TOPMED_MAX_PARALLEL_JOBS
    label "bigmem"

    input:
        path topmed_api_token_file
        path genotypes

    output:
        path "*.txt", emit: jobs_files

    script:
        genotypes_prefix = get_prefix(genotypes[0])
        """
        ${get_julia_cmd(task.cpus)} impute \
            genomicc \
            ${topmed_api_token_file} \
            --password ${params.TOPMED_ENCRYPTION_PASSWORD} \
            --max-concurrent-submissions ${params.TOPMED_MAX_PARALLEL_JOBS} \
            --refresh-rate ${params.TOPMED_REFRESH_RATE} \
            --r2 ${params.IMPUTATION_R2_FILTER} \
            --output-prefix genomicc
        """
}

process GetTOPMedDownloadList {
    input:
        path topmed_api_token_file
        val job_id

    output:
        tuple val(job_id), path("*.txt"), emit: info_files
        tuple val(job_id), path("*.zip"), emit: zip_files
        tuple val(job_id), path("*.md5", arity: 1), emit: md5_file
        tuple val(job_id), path("*.html", arity: 1), emit: report_file

    script:
        """
        ${get_julia_cmd(task.cpus)} get-topmed-download-list \
            ${job_id} \
            ${topmed_api_token_file} \
            --refresh-rate ${params.TOPMED_REFRESH_RATE}
        """
}

process DownloadTOPMedZipFile {
    input:
        tuple val(job_id), path(md5_file), path(zip_file_info)
        path topmed_api_token_file

    output:
        path output_file

    script:
        def output_file_parts = zip_file_info.getName().tokenize(".")
        output_file_parts.remove(1)
        output_file = output_file_parts.join(".")
        """
        ${get_julia_cmd(task.cpus)} download-topmed-file \
            ${job_id} \
            ${topmed_api_token_file} \
            ${zip_file_info} \
            --md5-file ${md5_file} \
            --refresh-rate ${params.TOPMED_REFRESH_RATE}
        """
}

process UnzipTOPMedFile {
    input:
        path(zip_file)

    output:
        path("*.dose.vcf.gz"), emit: vcf_file
        path("*.empiricalDose.vcf.gz"), emit: vcf_index_file
        path("*.info.txt"), emit: info_file
        path("*.html"), emit: report_file

    script:
        """
        unzip $zip_file
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

workflow Imputation {
    topmed_api_token = file(params.TOPMED_TOKEN_FILE)
    bed_genotypes = Channel.fromPath("${params.GENOTYPES_PREFIX}.{bed,bim,fam}").collect()

    if (params.TOPMED_JOBS_LIST == "NO_TOPMED_JOBS") {
        split_files = WriteImputationSplitLists(bed_genotypes)
        chrs_samples_split_files = split_files.chromosomes.splitText(){x -> x[0..-2]}
            .combine(split_files.samples.flatten())
        vcf_splits = MakeVCFSplit(
            bed_genotypes,
            chrs_samples_split_files
        )
        jobs_files = TOPMedImputation(topmed_api_token, vcf_splits.collect())
        job_ids = jobs_files.splitText().map {it -> it[0..-2]}
    }
    else {
        job_ids = Channel.fromList(params.TOPMED_JOBS_LIST)
    }

    files_to_download = GetTOPMedDownloadList(topmed_api_token, job_ids)

    zip_files_infos = files_to_download.zip_files.transpose()
    md5_files = files_to_download.md5_file.transpose()
    md5_to_zip_files_infos = md5_files.combine(zip_files_infos, by: 0)
    zip_files = DownloadTOPMedZipFile(md5_to_zip_files_infos, topmed_api_token)


    // // Make BGEN Output
    // // Idea 1: Merge vcf and then convert
    // all_dose_vcfs_indices = IndexVCF(all_dose_vcfs)
    // dose_vcfs_and_indices_by_chr = all_dose_vcfs.concat(all_dose_vcfs_indices)
    //     .map{ it -> [it.getName().tokenize(".")[1], it]}
    //     .groupTuple()
    //     .view()
    // chr_vcfs = MergeVCFsByChr(dose_vcfs_and_indices_by_chr)
    // pgen_files = VCFToPGEN(chr_vcfs)



}
