process CombineCovariates {
    publishDir "${params.PUBLISH_DIR}/covariates", mode: 'copy'

    input:
        path covariates_file
        path ancestry_file
        tuple path(eigenvalues), path(eigenvectors)

    output:
        path("covariates.merged.csv")

    script:
        """
        ${params.JULIA_CMD} combine-covariates \
            ${covariates_file} \
            ${ancestry_file} \
            ${eigenvectors} \
            --output=covariates.merged.csv
        """
}