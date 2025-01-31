include { get_prefix} from './utils.nf'

process QCRawGenotypes{
    publishDir "results/array-genotypes/qced", mode: 'symlink'

    input:
        tuple path(bed_file), path(bim_file), path(fam_file)
        path variants_to_flip

    output:
        tuple path("${output_prefix}.bed"), path("${output_prefix}.bim"), path("${output_prefix}.fam"), emit: genotypes 
        tuple path("${output_prefix}.filtered_samples.csv"), path("${output_prefix}.filtered_variants.csv"), emit: reports

    script:
        input_prefix = get_prefix(bed_file)
        output_prefix = "${input_prefix}.qced"
        """
        /opt/miniforge3/bin/mamba run -n plink_env plink --bfile ${input_prefix} \
            --noweb \
            --geno ${params.QC_GENOTYPE_MISSING_RATE} \
            --mind ${params.QC_INDIVIDUAL_MISSING_RATE} \
            --hwe ${params.QC_HWE} \
            --flip ${variants_to_flip} \
            --make-bed \
            --out ${input_prefix}.qced
        julia --project=/opt/sequential-gwas/ /opt/sequential-gwas/bin/seq-gwas.jl \
            report-qc-effect \
            ${input_prefix} ${output_prefix}
        """
}

process QCMergedGenotypes {
    publishDir "results/array-genotypes/merged_qced", mode: 'symlink'

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
            --noweb \
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