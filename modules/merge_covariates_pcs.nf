include { get_julia_cmd } from './utils.nf'

process MergeCovariatesPCs {
    label "bigmem"
    publishDir "${params.PUBLISH_DIR}/gwas/${group}", mode: 'symlink'

    input:
        tuple val(group), path(covariates), path(eigenvec), path(eigenval)

    output:
        tuple val(group), path("${output_file}")

    script:
        output_file = "${group}.covariates_pcs.csv"
        """
        ${get_julia_cmd(task.cpus)} merge-covariates-pcs \
            ${covariates} \
            ${eigenvec} \
            --output=${output_file}
        """
}