include { get_prefix } from './utils.nf'

process MakePopFile {
    publishDir "${params.ANCESTRY_PUBLISH_DIR}/popfile", mode: 'symlink'

    input:
        tuple path(bed_file), path(bim_file), path(fam_file)
        path(pedigree_file)

    output:
        path("${input_prefix}.pop")

    script:
        input_prefix = get_prefix(bed_file)
        """
        julia --project=/opt/sequential-gwas/ /opt/sequential-gwas/bin/seq-gwas.jl \
            make-pop-file \
            ${fam_file} \
            ${pedigree_file} \
            --output ${input_prefix}.pop
        """

}