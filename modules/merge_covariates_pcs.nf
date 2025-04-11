include { get_julia_cmd } from './utils.nf'

process MergeCovariatesPCs {
    label "bigmem"
    publishDir "${params.PUBLISH_DIR}/gwas/covariates", mode: 'symlink'

    input:
        path covariates
        path eigenvecs

    output:
        path("${output_file}")

    script:
        output_file = "covariates_and_pcs.csv"
        """
        ${get_julia_cmd(task.cpus)} merge-covariates-pcs \
            ${covariates} \
            chr \
            --output=${output_file}
        """
}