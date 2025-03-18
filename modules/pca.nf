include { get_prefix } from './utils.nf'

process GroupPCA {
    label "bigmem"
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
            --threads ${task.cpus} \
            --memory ${task.memory.toMega().toString()} \
            --bfile ${input_prefix} \
            --pca ${params.N_PCS} approx \
            --out ${group}
        """
}