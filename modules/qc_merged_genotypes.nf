include { get_prefix; get_julia_cmd } from './utils.nf'

process QCMergedGenotypes {
    label "multithreaded"
    publishDir "${params.MERGED_PUBLISH_DIR}/qced", mode: 'symlink'

    input:
        tuple path(bed_file), path(bim_file), path(fam_file)
        path unrelated_samples

    output:
        tuple path("${output_prefix}.bed"), path("${output_prefix}.bim"), path("${output_prefix}.fam"), emit: genotypes 
        path("${output_prefix}.log"), emit: log

    script:
        input_prefix = get_prefix(bed_file)
        output_prefix = "${input_prefix}.qced"
        """
        plink2 \
            --threads ${task.cpus} \
            --memory ${task.memory.toMega().toString()} \
            --bfile ${input_prefix} \
            --geno ${params.QC_GENOTYPE_MISSING_RATE} \
            --mind ${params.QC_INDIVIDUAL_MISSING_RATE} \
            --hwe ${params.QC_HWE_P} ${params.QC_HWE_K} \
            --output-chr chr26 \
            --keep ${unrelated_samples} \
            --make-bed \
            --out ${output_prefix}
        """
}