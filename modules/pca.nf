include { get_prefix } from './utils.nf'

process LocoGroupPCA {
    label "bigmem"
    label "multithreaded"
    publishDir "${params.PUBLISH_DIR}/gwas/${group}/pca", mode: 'symlink'

    input:
        tuple val(chr), val(group), path(bed_file), path(bim_file), path(fam_file)
    
    output:
        tuple val(chr), val(group), path("${output_prefix}.eigenval"), path("${output_prefix}.eigenvec")
    
    script:
        input_prefix = get_prefix(bed_file)
        approx = params.PCA_APPROX == true ? " approx" : ""
        output_prefix = "${chr}.${group}"
        """
        plink2 \
            --threads ${task.cpus} \
            --memory ${task.memory.toMega().toString()} \
            --bfile ${input_prefix} \
            --not-chr ${chr.replace('chr', '')} \
            --pca ${params.N_PCS}${approx} \
            --out ${output_prefix}
        """
}