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

process BasicQC{
    input:
        tuple path(bed_file), path(bim_file), path(fam_file)

    output:
        tuple path(bed_file), path(bim_file), path(fam_file)

    script:
    bed_prefix = get_prefix(bed_file)
        """
        /opt/miniforge3/bin/mamba run -n plink2_env plink2 --bfile ${bed_prefix} \
            --geno ${params.QC_GENOTYPE_MISSING_RATE} \
            --mind ${params.QC_INDIVIDUAL_MISSING_RATE} \ 
            --maf ${params.QC_MAF} \
            --hwe ${params.QC_HWE} \
            --make-bed 
            --out ${bed_prefix}.qced
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
    chain_file = Channel.fromPath(params.GRC37_TO_GRC38_CHAIN_FILE, checkIfExists: true)
    r8_array = Channel.fromFilePairs(params.R8_GENOTYPES, checkIfExists: true)
    before_2024_array = Channel.fromFilePairs(params.BEFORE_2024_GENOTYPES, checkIfExists: true)
    since_2024_array = Channel.fromFilePairs(params.SINCE_2024_GENOTYPES, checkIfExists: true)

    grc38_genotypes = PedToBed(since_2024_array)
    grc37_genotypes = PedToBed(r8_array.concat(before_2024_array))
    qced_grc38_genotypes = BasicQC(grc38_genotypes)
    qced_grc37_genotypes = BasicQC(grc37_genotypes)
    LiftOver(qced_grc37_genotypes, chain_file)
}