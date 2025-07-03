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
        path("ref_alleles.txt"), emit: ref_allele

    script:
        input_prefix = get_prefix(bed_file)
        output_prefix = "${input_prefix}.qced"
        """
        awk 'BEGIN {OFS="\\t"} {split(\$2, arr, ":"); print \$2, arr[3]}' ${bim_file} > ref_alleles.txt

        plink2 \
            --threads ${task.cpus} \
            --memory ${task.memory.toMega().toString()} \
            --bfile ${input_prefix} \
            --geno ${params.QC_GENOTYPE_MISSING_RATE} \
            --ref-allele ref_alleles.txt \
            --output-chr chr26 \
            --keep ${unrelated_samples} \
            --make-bed \
            --out ${output_prefix}

        awk '{ \$1 = \$2; print }' ${output_prefix}.fam > temp.fam
        mv temp.fam ${output_prefix}.fam
        """
}