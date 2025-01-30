include { get_prefix } from '../modules/utils.nf'

process PedToBed {
    publishDir "results/array-genotypes/bed", mode: 'symlink'

    input:
        tuple path(map_file), path(ped_file)

    output:
        tuple path("${input_prefix}.bed"), path("${input_prefix}.bim"), path("${input_prefix}.fam")

    script:
        input_prefix = get_prefix(ped_file)
        """
        /opt/miniforge3/bin/mamba run -n plink_env plink --file ${input_prefix} --make-bed --out ${input_prefix}
        """
}

process QCRawGenotypes{
    publishDir "results/array-genotypes/qced", mode: 'symlink'

    input:
        tuple path(bed_file), path(bim_file), path(fam_file)
        path variants_to_flip

    output:
        tuple path("${input_prefix}.qced.bed"), path("${input_prefix}.qced.bim"), path("${input_prefix}.qced.fam")

    script:
        input_prefix = get_prefix(bed_file)
        """
        /opt/miniforge3/bin/mamba run -n plink_env plink --bfile ${input_prefix} \
            --geno ${params.QC_GENOTYPE_MISSING_RATE} \
            --mind ${params.QC_INDIVIDUAL_MISSING_RATE} \
            --hwe ${params.QC_HWE} \
            --flip ${variants_to_flip} \
            --make-bed \
            --out ${input_prefix}.qced
        """
}

process LiftOver {
    publishDir "results/array-genotypes/lifted_over", mode: 'symlink'

    input:
        tuple path(map_file), path(ped_file)
        path chain_file
        
    output:
        tuple path("${input_prefix}.liftedOver.map"), path("${input_prefix}.liftedOver.ped"), emit: genotypes
        path("${input_prefix}.bed.unlifted"), optional: true

    script:
        input_prefix = get_prefix(map_file)
        """
        /opt/miniforge3/bin/mamba run -n liftover_env python /opt/sequential-gwas/bin/liftover.py \
            -m ${input_prefix}.map \
            -p ${input_prefix}.ped \
            -o ${input_prefix}.liftedOver \
            -c ${chain_file}
        """
}

def make_merge_list(merge_list, genotype_files) {
    def bed_files = genotype_files.findAll() { it.toString().endsWith(".bed") }
    println(bed_files)
    bed_files.each { it -> merge_list << "${get_prefix(it)}\n" }
}

process MergeGenotypes {
    publishDir "results/array-genotypes/merged", mode: 'symlink'

    input:
        path genotype_files
        path merge_list

    output:
        tuple path("genotypes.merged.bed"), path("genotypes.merged.bim"), path("genotypes.merged.fam")

    script:
        """
        /opt/miniforge3/bin/mamba run -n plink_env plink \
            --merge-list ${merge_list} \
            --make-bed \
            --out genotypes.merged
        """
}

workflow GenotypesQC {
    take: 
        grc37_genotypes
        grc38_genotypes
        variants_to_flip
        chain_file

    main:
        lifted_genotypes = LiftOver(grc37_genotypes, chain_file)
        genotypes_bed = PedToBed(lifted_genotypes.genotypes.concat(grc38_genotypes))
        qced_genotypes = QCRawGenotypes(genotypes_bed, variants_to_flip)
    
    emit:
        qced_genotypes
}

process QCMergedGenotypes {
    input:
        path genotypes

    output:
        tuple path("genotypes.merged.qced.bed"), path("genotypes.merged.qced.bim"), path("genotypes.merged.qced.fam")

    script:
        """
        /opt/miniforge3/bin/mamba run -n plink_env plink --bfile genotypes.merged \
            --geno ${params.QC_GENOTYPE_MISSING_RATE} \
            --mind ${params.QC_INDIVIDUAL_MISSING_RATE} \
            --hwe ${params.QC_HWE} \
            --make-bed \
            --out genotypes.merged.qced
        """
}

workflow CombineGenotypingArrays {
    // GRC37 Genotypes
    r8_array = Channel.fromPath("${params.R8_GENOTYPES}*", checkIfExists: true).collect()
    before_2024_array = Channel.fromPath("${params.BEFORE_2024_GENOTYPES}*", checkIfExists: true).collect()
    grc37_genotypes = r8_array.concat(before_2024_array)
    // GRC38 Genotypes
    grc38_genotypes = Channel.fromPath("${params.SINCE_2024_GENOTYPES}*", checkIfExists: true).collect()
    // LiftOver and QC Genotypes
    variants_to_flip = file(params.VARIANTS_TO_FLIP_GRC38, checkIfExists: true)
    chain_file = file(params.GRC37_TO_GRC38_CHAIN_FILE, checkIfExists: true)
    qced_genotypes = GenotypesQC(grc37_genotypes, grc38_genotypes, variants_to_flip, chain_file)
    // Merge Genotypes
    merge_list = qced_genotypes
        .map { it -> get_prefix(it[0].getName()) }
        .collectFile(name: "merge_list.txt", newLine: true)
    merged_genotypes = MergeGenotypes(qced_genotypes.collect(), merge_list)
    QCMergedGenotypes(merged_genotypes)
}