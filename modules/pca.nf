include { get_prefix } from './utils.nf'

process GroupPCA {
    label "multithreaded"
    publishDir "${params.PUBLISH_DIR}/gwas/${group}/pca", mode: 'symlink'

    input:
        tuple val(group), path(bed_file), path(bim_file), path(fam_file)
    
    output:
        tuple val(group), path("${group}.eigenvec"), path("${group}.eigenval")
    
    script:
        input_prefix = get_prefix(bed_file)
        """
        plink2 \
            --bfile ${input_prefix} \
            --pca ${params.N_PCS} \
            --out ${group}
        """
}