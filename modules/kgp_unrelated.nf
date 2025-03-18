include { get_prefix; get_julia_cmd } from './utils.nf'

process KeepKGPUnrelated {
    label "multithreaded"
    label "bigmem"
    publishDir "${params.KGP_PUBLISH_DIR}/merged_unrelated", mode: 'symlink'

    input:
        tuple path(bed_file), path(fam_file), path(bim_file)
        path pedigree_file

    output:
        tuple path("${output_prefix}.bed"), path("${output_prefix}.bim"), path("${output_prefix}.fam"), path("${output_prefix}.acount"), emit: genotypes
        path "kgp_unrelated_individuals.txt", emit: unrelated_individuals

    script:
        // Keep only unrelated individuals from the KGP dataset
        input_prefix = get_prefix(bed_file)
        output_prefix = "${input_prefix}.unrelated"
        """
        ${get_julia_cmd(task.cpus)} get-kgp-unrelated-individuals \
            ${pedigree_file} \
            --output kgp_unrelated_individuals.txt

        plink2 \
            --threads ${task.cpus} \
            --memory ${task.memory.toMega().toString()} \
            --bfile ${input_prefix} \
            --keep kgp_unrelated_individuals.txt \
            --output-chr chr26 \
            --freq counts \
            --make-bed \
            --out ${output_prefix}
        """

}