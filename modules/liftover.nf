include { get_prefix } from './utils.nf'

process LiftOver {
    publishDir "${params.ARRAY_GENOTYPES_PUBLISH_DIR}/lifted_over", mode: 'symlink'

    input:
        tuple val(id), path(ped_file), path(map_file)
        path chain_file
        
    output:
        tuple val(id), path("${output_prefix}.ped"), path("${output_prefix}.map"), emit: genotypes
        tuple val(id), path("${output_prefix}.bed.unlifted"), emit: unlifted

    script:
        input_prefix = get_prefix(map_file)
        output_prefix = "${input_prefix}.liftedOver"
        """
        /opt/miniforge3/bin/mamba run -n liftover_env python /opt/sequential-gwas/bin/liftover.py \
            -m ${input_prefix}.map \
            -p ${input_prefix}.ped \
            -o ${input_prefix}.liftedOver \
            -c ${chain_file}
        """
}