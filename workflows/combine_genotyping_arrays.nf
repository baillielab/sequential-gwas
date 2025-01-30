include { get_prefix } from '../modules/utils.nf'

process PedToBed {
    input:
        tuple val(genome_build), path(map_file), path(ped_file)

    output:
        tuple val(genome_build), path("${input_prefix}.bed"), path("${input_prefix}.bim"), path("${input_prefix}.fam")

    script:
        input_prefix = get_prefix(ped_file)
        """
        /opt/miniforge3/bin/mamba run -n plink_env plink --file ${input_prefix} --make-bed --out ${input_prefix}
        """
}

process BasicQC{
    input:
        tuple val(genome_build), path (variants_to_flip), path(bed_file), path(bim_file), path(fam_file)

    output:
        tuple val(genome_build), path("${input_prefix}.qced.bed"), path("${input_prefix}.qced.bim"), path("${input_prefix}.qced.fam")

    script:
        input_prefix = get_prefix(bed_file)
        """
        /opt/miniforge3/bin/mamba run -n plink_env plink --bfile ${input_prefix} \
            --geno ${params.QC_GENOTYPE_MISSING_RATE} \
            --mind ${params.QC_INDIVIDUAL_MISSING_RATE} \
            --maf ${params.QC_MAF} \
            --hwe ${params.QC_HWE} \
            --flip ${variants_to_flip} \
            --make-bed \
            --out ${input_prefix}.qced
        """
}

process LiftOver {
    publishDir "results/lifted_over_genotypes", mode: 'symlink'

    input:
        tuple val(genome_build), path(chain_file), path(map_file), path(ped_file)
        
    output:
        tuple val(genome_build), path("${input_prefix}.liftedOver.map"), path("${input_prefix}.liftedOver.ped"), emit: genotypes
        path("${input_prefix}.bed.unlifted"), emit: unlifted

    script:
        input_prefix = get_prefix(map_file)
        if (chain_file.getName() != "NO_LIFTOVER_NEEDED") { 
            """
            /opt/miniforge3/bin/mamba run -n liftover_env python /mnt/code/bin/liftover.py \
                -m ${input_prefix}.map \
                -p ${input_prefix}.ped \
                -o ${input_prefix}.liftedOver \
                -c ${chain_file}
            """
        }
        else {
            """
            ln -s ${input_prefix}.ped ${input_prefix}.liftedOver.ped
            ln -s ${input_prefix}.map ${input_prefix}.liftedOver.map
            """
        }
}

workflow GenotypesQC {
    take: 
        genotypes
        variants_to_flip
        chain_file

    main:
        lifted_genotypes = LiftOver(chain_file.combine(genotypes, by: 0))
        genotypes_bed = PedToBed(lifted_genotypes.genotypes)
        qced_genotypes = BasicQC(variants_to_flip.combine(genotypes_bed, by: 0))
    
    emit:
        qced_genotypes
}

workflow CombineGenotypingArrays {
    // GRC37 Genotypes
    r8_array = Channel.fromPath("${params.R8_GENOTYPES}*", checkIfExists: true).collect()
    before_2024_array = Channel.fromPath("${params.BEFORE_2024_GENOTYPES}*", checkIfExists: true).collect()
    grc37_genotypes = r8_array.concat(before_2024_array)
    // GRC38 Genotypes
    grc38_genotypes = Channel.fromPath("${params.SINCE_2024_GENOTYPES}*", checkIfExists: true).collect()
    // QC and LiftOver Genotypes
    genotypes = Channel.of("grch37", "grch37", "grch38")
        .merge(grc37_genotypes.concat(grc38_genotypes))
    variants_to_flip = Channel.of("grch37", "grch38")
        .merge(Channel.fromPath([params.VARIANTS_TO_FLIP_GRC37, params.VARIANTS_TO_FLIP_GRC38], checkIfExists: true))
    chains_files = Channel.of("grch37", "grch38")
        .merge(Channel.fromPath([params.GRC37_TO_GRC38_CHAIN_FILE, "${projectDir}/assets/NO_LIFTOVER_NEEDED"], checkIfExists: true))
    qced_liftedover_genotypes = GenotypesQC(genotypes, variants_to_flip, chains_files)

}