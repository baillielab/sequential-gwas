include { get_prefix } from '../modules/utils.nf'

process Admixture {
    input:
        path genotypes
        val k
    output:
        tuple path("${input_prefix}.${k}.P"), path("${input_prefix}.${k}.Q"), emit: pq_files
        path("${input_prefix}.${k}.log"), emit: log_file
    script:
        input_prefix = get_prefix(genotypes)
        """
        admixture ${input_prefix}.bed --cv genotypes.bed ${k} -j${task.cpu} -s 12345 | tee ${input_prefix}.${k}.admixturelog
        """
}

process FindNumberOfPopulations{
    input:
        file admixture_log_files
    output:
        val "${outputfile.baseName}"
    script:
        """
        grep -h CV *.admixturelog > admixture_cross_validation_results.txt
        """
}

workflow EstimateAncestry {
    take:
        genotypes
        reference_panel
    main:
        ancestry = EstimateAncestry(genotypes, reference_panel)
    emit:
        ancestry
}