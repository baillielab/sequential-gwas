include { get_prefix } from './utils.nf'

process FilterHighLoadingsVariants {
    label "multithreaded"
    publishDir "${params.PUBLISH_DIR}", mode: 'copy'

    input:
        tuple path(bed_file), path(bim_file), path(fam_file)
        path(high_loadings_variants)

    output:
        tuple path("${output_prefix}.bed"), path("${output_prefix}.bim"), path("${output_prefix}.fam")

    script:
        input_prefix = get_prefix(bed_file)
        output_prefix = "genotypes.aggregated.qced.final"
        if (params.FILTER_HIGH_LOADINGS_VARIANTS == true)
            """
            plink2 \
                --threads ${task.cpus} \
                --memory ${task.memory.toMega().toString()} \
                --bfile ${input_prefix} \
                --output-chr chr26 \
                --exclude ${high_loadings_variants} \
                --make-bed \
                --out ${output_prefix}
            """
        else // This copying is not ideal but makes it easy to publish to the same directory
            """cp ${input_prefix}.bed ${output_prefix}.bed
               cp ${input_prefix}.bim ${output_prefix}.bim
               cp ${input_prefix}.fam ${output_prefix}.fam
            """
        

}