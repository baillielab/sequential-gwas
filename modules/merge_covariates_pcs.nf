process MergeCovariatesPCs {
    publishDir "${params.PUBLISH_DIR}/gwas/${group}", mode: 'symlink'

    input:
        tuple val(group), path(covariates), path(eigenvec), path(eigenval)

    output:
        path "${output_file}"

    script:
        output_file = "${group}.covariates_pcs.csv"
        """
        ${params.JULIA_CMD} merge-covariates-pcs \
            ${covariates} \
            ${eigenvec} \
            --output=${output_file}
        """
}