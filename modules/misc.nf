include { get_prefix} from './utils.nf'

process PedToBed {
    publishDir "results/array-genotypes/bed", mode: 'symlink'

    input:
        tuple path(map_file), path(ped_file)

    output:
        tuple path("${input_prefix}.bed"), path("${input_prefix}.bim"), path("${input_prefix}.fam")

    script:
        input_prefix = get_prefix(ped_file)
        """
        /opt/miniforge3/bin/mamba run -n plink_env plink \
            --noweb \
            --alleleACGT \
            --file ${input_prefix} \
            --make-bed \
            --out ${input_prefix}
        """
}

process LiftOver {
    publishDir "results/array-genotypes/lifted_over", mode: 'symlink'

    input:
        tuple path(map_file), path(ped_file)
        path chain_file
        
    output:
        tuple path("${input_prefix}.liftedOver.map"), path("${input_prefix}.liftedOver.ped"), emit: genotypes
        path("${input_prefix}.liftedOver.bed.unlifted"), optional: true

    script:
        input_prefix = get_prefix(map_file)
        """
        /opt/miniforge3/bin/mamba run -n liftover_env python /opt/sequential-gwas/bin/liftover.py \
            -m ${input_prefix}.map \
            -p ${input_prefix}.ped \
            -o ${input_prefix}.liftedOver \
            -c ${chain_file}
        """
}

process MergeGenotypes {
    publishDir "results/array-genotypes/merged", mode: 'symlink'

    input:
        path genotype_files
        path merge_list

    output:
        tuple path("genotypes.merged.bed"), path("genotypes.merged.bim"), path("genotypes.merged.fam")

    script:
        """
        /opt/miniforge3/bin/mamba run -n plink_env plink \
            --noweb \
            --merge-list ${merge_list} \
            --make-bed \
            --out genotypes.merged
        """
}