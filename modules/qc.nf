include { get_prefix} from './utils.nf'

process QCRawGenotypes{
    publishDir "${params.ARRAY_GENOTYPES_PUBLISH_DIR}/qced", mode: 'symlink'

    input:
        tuple val(id), path(bed_file), path(bim_file), path(fam_file)
        path variants_to_flip

    output:
        tuple val(id), path("${output_prefix}.bed"), path("${output_prefix}.bim"), path("${output_prefix}.fam"), emit: genotypes 
        tuple path("${output_prefix}.filtered_samples.csv"), path("${output_prefix}.filtered_variants.csv"), emit: reports

    script:
        input_prefix = get_prefix(bed_file)
        output_prefix = "${input_prefix}.qced"
        """
        /opt/miniforge3/bin/mamba run -n plink2_env plink2 --bfile ${input_prefix} \
            --geno ${params.QC_GENOTYPE_MISSING_RATE} \
            --mind ${params.QC_INDIVIDUAL_MISSING_RATE} \
            --hwe ${params.QC_HWE} \
            --set-all-var-ids @:# \
            --output-chr chr26 \
            --make-bed \
            --out ${input_prefix}.qced
        julia --project=/opt/sequential-gwas/ /opt/sequential-gwas/bin/seq-gwas.jl \
            report-qc-effect \
            ${input_prefix} ${output_prefix}
        """
}

process QCMergedGenotypes {
    publishDir "${params.ARRAY_GENOTYPES_PUBLISH_DIR}/merged_qced", mode: 'symlink'

    input:
        tuple path(bed_file), path(bim_file), path(fam_file)

    output:
        tuple path("${output_prefix}.bed"), path("${output_prefix}.bim"), path("${output_prefix}.fam"), emit: genotypes 
        tuple path("${output_prefix}.filtered_samples.csv"), path("${output_prefix}.filtered_variants.csv"), emit: reports

    script:
        input_prefix = get_prefix(bed_file)
        output_prefix = "${input_prefix}.qced"
        """
        /opt/miniforge3/bin/mamba run -n plink_env plink --bfile ${input_prefix} \
            --geno ${params.QC_GENOTYPE_MISSING_RATE} \
            --mind ${params.QC_INDIVIDUAL_MISSING_RATE} \
            --hwe ${params.QC_HWE} \
            --make-bed \
            --out ${output_prefix}
        julia --project=/opt/sequential-gwas/ /opt/sequential-gwas/bin/seq-gwas.jl \
            report-qc-effect \
            ${input_prefix} ${output_prefix}
        """
}