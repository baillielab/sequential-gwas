process KingRelationshipInference {
    label "multithreaded"
    publishDir "${params.MERGED_PUBLISH_DIR}/king_relatedness", mode: 'symlink'

    input:
        tuple path(bed_file), path(bim_file), path(fam_file)

    output:
        path("kingunrelated.txt")

    script:
        """
        sed 's/^chr//' ${bim_file} > no_chr.bim
        king --cpus ${task.cpus} -b ${bed_file} --bim no_chr.bim --fam ${fam_file} --unrelated --degree 2
        """
}