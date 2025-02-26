include { get_prefix } from './utils.nf'

process QCMergedGenotypes {
    publishDir "${params.MERGED_PUBLISH_DIR}/qced", mode: 'symlink'

    input:
        tuple path(bed_file), path(bim_file), path(fam_file)
        path unrelated_samples

    output:
        tuple path("${output_prefix}.bed"), path("${output_prefix}.bim"), path("${output_prefix}.fam"), emit: genotypes 
        tuple path("${output_prefix}.filtered_samples.csv"), path("${output_prefix}.filtered_variants.csv"), emit: reports

    script:
        input_prefix = get_prefix(bed_file)
        output_prefix = "${input_prefix}.qced"
        """
        plink --bfile ${input_prefix} \
            --geno ${params.QC_GENOTYPE_MISSING_RATE} \
            --mind ${params.QC_INDIVIDUAL_MISSING_RATE} \
            --hwe ${params.QC_HWE_P} \
            --output-chr chr26 \
            --keep ${unrelated_samples} \
            --make-bed \
            --out ${output_prefix}
        
        ${params.JULIA_CMD} report-qc-effect \
            ${input_prefix} ${output_prefix}
        """
}