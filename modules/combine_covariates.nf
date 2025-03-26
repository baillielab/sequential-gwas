include { get_julia_cmd } from './utils.nf'

process CombineCovariates {
    label "bigmem"
    publishDir "${params.PUBLISH_DIR}", mode: 'copy'

    input:
        path ancestry_file
        tuple path(eigenvalues), path(eigenvectors)
        path wgs_samples
        path release_r8_fam
        path release_2021_2023_fam
        path release_2024_now_fam

    output:
        path("covariates.merged.csv")

    script:
        """
        ${get_julia_cmd(task.cpus)} combine-covariates \
            ${ancestry_file} \
            ${eigenvectors} \
            ${wgs_samples} \
            ${release_r8_fam} \
            ${release_2021_2023_fam} \
            ${release_2024_now_fam} \
            --output=covariates.merged.csv
        """
}