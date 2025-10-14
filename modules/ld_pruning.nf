include { get_prefix } from './utils.nf'

process LDPruning {
    label "multithreaded"
    
    input:
        tuple path(bed_file), path(bim_file), path(fam_file)
        path high_ld_regions

    output:
        tuple path("${output_prefix}.bed"), path("${output_prefix}.bim"), path("${output_prefix}.fam")

    script:
        input_prefix = get_prefix(bed_file)
        output_prefix = "${input_prefix}.ldpruned"
        """
        plink2 \
            --threads ${task.cpus} \
            --memory ${task.memory.toMega().toString()} \
            --bfile ${input_prefix} \
            --indep-pairwise ${params.IP_VALUES}
        
        plink2 \
            --threads ${task.cpus} \
            --memory ${task.memory.toMega().toString()} \
            --bfile ${input_prefix} \
            --extract plink2.prune.in \
            --maf ${params.PCA_MAF} \
            --make-bed \
            --exclude range ${high_ld_regions} \
            --out ${output_prefix} \
        """
}