include { get_prefix } from './utils.nf'

process LiftOver {
    publishDir "${params.ARRAY_GENOTYPES_PUBLISH_DIR}/lifted_over", mode: 'symlink'

    input:
        tuple val(id), path(ped_file), path(map_file)
        path chain_file
        
    output:
        tuple val(id), path("${output_prefix}.ped"), path("${output_prefix}.map"), emit: genotypes
        path("${input_prefix}.liftover_report.csv"), emit: report

    script:
        input_prefix = get_prefix(map_file)
        output_prefix = "${input_prefix}.liftedOver"
        """
        /opt/miniforge3/bin/mamba run -n liftover_env python /opt/sequential-gwas/bin/liftover.py \
            -m ${input_prefix}.map \
            -p ${input_prefix}.ped \
            -o ${input_prefix}.liftedOver \
            -c ${chain_file}
        julia --project=/opt/sequential-gwas/ /opt/sequential-gwas/bin/seq-gwas.jl \
            report-liftover-effect \
            ${input_prefix} ${output_prefix}
        """
}