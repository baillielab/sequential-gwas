include { get_prefix} from './utils.nf'

process PedToBed {
    publishDir "${params.ARRAY_GENOTYPES_PUBLISH_DIR}/bed", mode: 'symlink'

    input:
        tuple val(id), path(ped_file), path(map_file)

    output:
        tuple val(id), path("${input_prefix}.bed"), path("${input_prefix}.bim"), path("${input_prefix}.fam")

    script:
        input_prefix = get_prefix(ped_file)
        """
        plink \
            --alleleACGT \
            --file ${input_prefix} \
            --biallelic-only strict \
            --make-bed \
            --out ${input_prefix}
        """
}