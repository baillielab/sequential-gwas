include { get_prefix } from './utils.nf'

process EstimateAncestry {
    publishDir "${params.ANCESTRY_PUBLISH_DIR}/ancestry", mode: 'symlink'

    input:
        tuple path(bed_file), path(bim_file), path(fam_file)
        path pedigree

    output:
        path("${input_prefix}.*.{P,Q}"), emit: pq_files
        path("${output}"), emit: ancestry

    script:
        input_prefix = get_prefix(bed_file)
        output = "${input_prefix}.ancestry.csv"
        """
        seq-gwas estimate-ancestry \
            ${input_prefix} \
            ${pedigree} \
            --output=${output}
        """
}