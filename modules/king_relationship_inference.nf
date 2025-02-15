process KingRelationshipInference {
    input:
        tuple path(bed_file), path(bim_file), path(fam_file)

    output:
        path("kingunrelated.txt")

    script:
        """
        sed 's/^chr//' ${bim_file} > no_chr.bim
        king -b ${bed_file} --bim no_chr.bim --fam ${fam_file} --unrelated --degree 2
        """
}