include { get_prefix } from './utils.nf'

process PlinkPCA {
    publishDir "${params.PCA_PUBLISH_DIR}/pca", mode: 'symlink'

    input:
        tuple path(bed_file), path(bim_file), path(fam_file)
    
    output:
        tuple path("${input_prefix}.eigenvec"), path("${input_prefix}.eigenval")
    
    script:
        input_prefix = get_prefix(bed_file)
        """
        /opt/miniforge3/bin/mamba run -n plink2_env plink2 \
            --bfile ${input_prefix} \
            --pca ${params.N_PCS} \
            --out ${input_prefix}
        """
}