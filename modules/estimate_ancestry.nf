include { get_prefix; get_julia_cmd } from './utils.nf'

process EstimateAncestry {
    label "hyperthreaded"

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
        ${get_julia_cmd(task.cpus)} estimate-ancestry \
            ${input_prefix} \
            ${pedigree} \
            --output=${output} \
            --threshold=${params.ANCESTRY_THRESHOLD}
        """
}