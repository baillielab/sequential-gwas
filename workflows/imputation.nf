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
        path("*.dose.vcf.gz"), emit: dose
        path("*.empiricalDose.vcf.gz"), emit: empirical_dose
        path("*.info.gz"), emit: info

    script:
        jobname = zip_file.getName().tokenize(".")[0]
        """
        unzip -P ${params.TOPMED_ENCRYPTION_PASSWORD} ${zip_file} -d temp_extract

        for f in temp_extract/*; do
            mv "\$f" "./${jobname}.\$(basename "\$f")"
        done
        """
}

process MergeVCFsByChr {
    input:
        tuple val(chr), path(vcf_files)

    output:
        path("${output}")

    script:
        output = "${chr}.vcf.gz"
        sorted_vcf_files_string = vcf_files
            .findAll { x -> x.getName().endsWith("vcf.gz") }
            .sort{ x -> x.getName().tokenize("_")[1].toInteger() }
            .join("\n")
        """
        echo "${sorted_vcf_files_string}" > merge_list.txt

        bcftools merge --threads ${task.cpus} -o ${output} -O z -l merge_list.txt
        """
}

process IndexVCF {
    input:
        path vcf_file

    output:
        path "${vcf_file}.tbi"

    script:
        """
        bcftools index --tbi --threads ${task.cpus} ${vcf_file}
        """
}

process QCMergedImputedFile {
    input:
        path vcf_file

    output:
        path "${output}"

    script:
        output = "${get_prefix(get_prefix(vcf_file))}.qced.vcf.gz"
        """
        bcftools view -m2 -e '( R2 < ${params.IMPUTATION_R2_FILTER})' --threads ${task.cpus} -O z -o ${output} ${vcf_file}
        """
}

process VCFToPGEN {
    label "bigmem"
    publishDir "${params.PUBLISH_DIR}", mode: "copy"

    input:
        path vcf_file

    output:
        tuple path("${output_prefix}.pgen"), path("${output_prefix}.pvar"), path("${output_prefix}.psam")

    script:
        output_prefix = get_prefix(get_prefix(vcf_file))
        """
        plink2 --vcf ${vcf_file} --make-pgen --threads ${task.cpus} --out ${output_prefix}

        awk -F'\\t' 'BEGIN{OFS="\\t"} {sub(/.*@/, "", \$1); print}' ${output_prefix}.psam > ${output_prefix}.psam.tmp
        mv ${output_prefix}.psam.tmp ${output_prefix}.psam
        """
}


workflow Imputation {
    topmed_api_token = file(params.TOPMED_TOKEN_FILE)
    bed_genotypes = Channel.fromPath("${params.GENOTYPES_PREFIX}.{bed,bim,fam}").collect()
    // Send for Imputation or retrieve jobs list
    if (params.TOPMED_JOBS_LIST == "NO_TOPMED_JOBS") {
        split_files = WriteImputationSplitLists(bed_genotypes)
        chrs_samples_split_files = split_files.chromosomes.splitText(){x -> x[0..-2]}
            .combine(split_files.samples.flatten())
        vcf_splits = MakeVCFSplit(
            bed_genotypes,
            chrs_samples_split_files
        )
        jobs_files = TOPMedImputation(topmed_api_token, vcf_splits.collect())
        job_ids = jobs_files
            .flatten()
            .splitText()
            .map { it -> it.trim() }
    }
    else {
        job_ids = Channel.fromList(params.TOPMED_JOBS_LIST)
    }
    // Download TOPMed files
    files_to_download = GetTOPMedDownloadList(topmed_api_token, job_ids)
    zip_files_infos = files_to_download.zip_files.transpose()
    md5_files = files_to_download.md5_file.transpose()
    md5_to_zip_files_infos = md5_files.combine(zip_files_infos, by: 0)
    zip_files = DownloadTOPMedZipFile(md5_to_zip_files_infos, topmed_api_token)
    // Unzip TOPMed files
    unziped_files = UnzipTOPMedFile(zip_files)
    // Merge VCFs by chromosome
    imputed_files = unziped_files.dose.flatten()
    indices = IndexVCF(imputed_files)
    imputed_files_and_indices_by_chr = imputed_files
        .concat(indices)
        .map { it -> [it.getName().tokenize(".")[1], it]}
        .groupTuple()
    chr_vcfs = MergeVCFsByChr(imputed_files_and_indices_by_chr)
    // QC merged VCFs
    qced_vcfs_chrs = QCMergedImputedFile(chr_vcfs)
    // Convert VCFs to PGEN
    VCFToPGEN(qced_vcfs_chrs)
}
