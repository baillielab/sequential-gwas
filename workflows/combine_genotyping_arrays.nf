params.GRC37_TO_GRC38_CHAIN_FILE = "${projectDir}/assets/hg19ToHg38.over.chain"

include { get_prefix } from './modules/utils.nf'

process PedToBed {
    input:
        tuple path(ped_file), path(map_file)

    output:
        tuple path("${ped_prefix}.bed"), path("${ped_prefix}.bim"), path("${ped_prefix}.fam")

    script:
        ped_prefix = get_prefix(ped_file)
        """
        /opt/miniforge3/bin/mamba run -n plink2_env plink2 --file ${ped_prefix} --make-bed
        """
}


process LiftOver {
    publishDir "results/lifted_over_genotypes", mode: 'symlink'
    input:
        path genotypes
        path chain_file

    output:
        tuple path("genotypesLiftedOver.bed"), path("genotypesLiftedOver.bim"), path("genotypesLiftedOver.fam")

    script:
        """
        /opt/miniforge3/bin/mamba run -n liftover_env plink2 --bfile genotypes --lift ${chain_file} --make-bed --out genotypesLiftedOver
        """
}


workflow CombineGenotypingArrays {
    chain_file = Channel.fromPath(params.CHAIN_FILE)
    r8_array_ch = Channel.fromFilePairs(params.R8_GENOTYPES)
    before_2024_array_ch = Channel.fromFilePairs(params.BEFORE_2024_GENOTYPES)
    since_2024_array_ch = Channel.fromFilePairs(params.SINCE_2024_GENOTYPES)

    array_ch = r8_array_ch.concat(before_2024_array_ch).concat(since_2024_array_ch)
    PedToBed(all_ped_files)
    // LiftOver()
}