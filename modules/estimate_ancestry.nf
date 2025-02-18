include { get_prefix } from './utils.nf'

process EstimateAncestry {
    input:
        path genotypes
        path popfile

    output:
        tuple path("${input_prefix}.5.P"), path("${input_prefix}.5.Q"), emit: pq_files
        path("${input_prefix}.5.admixturelog"), emit: log_file

    script:
        input_prefix = get_prefix(genotypes)
        """
        admixture ${input_prefix}.bed --supervised -j${task.cpu} -s 123 | tee ${input_prefix}.5.admixturelog
        """
}