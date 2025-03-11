include { get_prefix } from './utils.nf'

process LiftOver {
    label "multithreaded"
    publishDir "${params.ARRAY_GENOTYPES_PUBLISH_DIR}/lifted_over", mode: 'symlink'

    input:
        tuple val(id), path(ped_file), path(map_file)
        path chain_file
        
    output:
        tuple val(id), path("${output_prefix}.bed"), path("${output_prefix}.bim"), path("${output_prefix}.fam"), emit: genotypes
        tuple val(id), path("${input_prefix}.liftover_temp.unmapped"), emit: unlifted

    script:
        input_prefix = get_prefix(map_file)
        output_prefix = "${input_prefix}.liftedOver"
        """
        /opt/miniforge3/bin/mamba run -n liftover_env python /opt/sequential-gwas/bin/liftover.py \
            -m ${map_file} \
            -o ${input_prefix}.liftover_temp \
            -c ${chain_file}
        plink \
            --ped ${input_prefix}.ped \
            --map ${input_prefix}.liftover_temp.map \
            --biallelic-only strict \
            --allow-extra-chr \
            --chr 1-22 \
            --alleleACGT \
            --exclude ${input_prefix}.liftover_temp.unmapped \
            --make-bed \
            --out ${output_prefix}
        """
}