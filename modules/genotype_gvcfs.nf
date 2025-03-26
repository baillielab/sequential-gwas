include { get_julia_cmd } from './utils.nf'

process GenotypeGVCFs{
    label "multithreaded"
    publishDir "${params.WGS_PUBLISH_DIR}/genotyped", mode: 'symlink'

    input:
        tuple val(prefix), path(gvcf)
        path shared_variants_gatk
        path shared_variants_plink
        path reference_genome
        path reference_genome_index
        path reference_genome_dict

    output:
        tuple path("${output_prefix}.bed"), path("${output_prefix}.bim"), path("${output_prefix}.fam")

    script:
        output_prefix = "${prefix}.shared"
        """
        ${get_julia_cmd(task.cpus)} genotype-gvcf \
            ${prefix}.gvcf.gz \
            ${shared_variants_plink} \
            ${shared_variants_gatk} \
            ${reference_genome} \
            --output ${output_prefix}
        """
}