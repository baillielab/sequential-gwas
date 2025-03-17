include { get_julia_cmd } from './utils.nf'

process CombineCovariates {
    publishDir "${params.PUBLISH_DIR}", mode: 'copy'

    input:
        path covariates_file
        path ancestry_file
        tuple path(eigenvalues), path(eigenvectors)

    output:
        path("covariates.merged.csv")

    script:
        """
        ${get_julia_cmd(task.cpus)} combine-covariates \
            ${covariates_file} \
            ${ancestry_file} \
            ${eigenvectors} \
            --output=covariates.merged.csv
        """
}