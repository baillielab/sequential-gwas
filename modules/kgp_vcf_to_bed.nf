include { get_prefix } from './utils.nf'

process KGPVCFToBed {
    publishDir "${params.KGP_PUBLISH_DIR}/bed", mode: 'symlink'

    input:
        tuple path(vcf_file), path(vcf_index)

    output:
        tuple path("${input_prefix}.bed"), path("${input_prefix}.bim"), path("${input_prefix}.fam")

    script:
        input_prefix = get_prefix(vcf_file)[0..-5]
        """
        plink2 \
            --output-chr chr26 \
            --vcf ${vcf_file} \
            --make-bed \
            --out ${input_prefix}
        """
}