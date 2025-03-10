include { get_prefix; get_julia_cmd } from './utils.nf'

process GenotypingArrayBasicQC {
    label "multithreaded"
    publishDir "${params.ARRAY_GENOTYPES_PUBLISH_DIR}/qced", mode: 'symlink'

    input:
        tuple val(id), path(bed_file), path(bim_file), path(fam_file)
        path variants_to_flip

    output:
        tuple val(id), path("${output_prefix}.bed"), path("${output_prefix}.bim"), path("${output_prefix}.fam"), path("${output_prefix}.acount"), emit: genotypes 
        tuple val(id), path("${output_prefix}.filtered_samples.csv"), path("${output_prefix}.filtered_variants.csv"), emit: reports

    script:
        input_prefix = get_prefix(bed_file)
        output_prefix = "${input_prefix}.qced"
        """
        plink2 --bfile ${input_prefix} \
            --geno ${params.QC_GENOTYPE_MISSING_RATE} \
            --mind ${params.QC_INDIVIDUAL_MISSING_RATE} \
            --hwe ${params.QC_HWE_P} ${params.QC_HWE_K} \
            --set-all-var-ids @:# \
            --output-chr chr26 \
            --freq counts \
            --make-bed \
            --out ${input_prefix}.qced
        
        ${get_julia_cmd(task.cpus)} report-qc-effect \
            ${input_prefix} ${output_prefix}
        """
}