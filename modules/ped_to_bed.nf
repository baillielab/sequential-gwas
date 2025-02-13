include { get_prefix} from './utils.nf'

process PedToBed {
    publishDir "${params.ARRAY_GENOTYPES_PUBLISH_DIR}/bed", mode: 'symlink'

    input:
        tuple path(map_file), path(ped_file)

    output:
        tuple path("${input_prefix}.bed"), path("${input_prefix}.bim"), path("${input_prefix}.fam")

    script:
        input_prefix = get_prefix(ped_file)
        """
        /opt/miniforge3/bin/mamba run -n plink_env plink \
            --alleleACGT \
            --file ${input_prefix} \
            --make-bed \
            --out ${input_prefix}
        """
}