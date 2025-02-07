include { get_prefix} from './utils.nf'

process MakeSharedVariantsList {
    publishDir "results/array-genotypes/qced_shared_variants", pattern: "variants_intersection.csv",  mode: 'symlink'
    publishDir "results/wgs/shared_variants", pattern: "variants_intersection.{bed,bim}", mode: 'symlink'

    input:
        path genotypes_files

    output:
        path("variants_intersection.csv"), emit: joint_variants_plink
        path("variants_intersection.bed"), emit: joint_variants_gatk
        path("variants_intersection.bim"), emit: joint_variants_bim

    script:
        """
        julia --project=/opt/sequential-gwas/ /opt/sequential-gwas/bin/seq-gwas.jl \
            write-variants-intersection
        """
}

process ExtractSharedVariantsFromPLINK {
    publishDir "results/array-genotypes/qced_shared_variants", mode: 'symlink'

    input:
        tuple path(bed_file), path(bim_file), path(fam_file)
        path joint_variants

    output:
        tuple path("${output_prefix}.bed"), path("${output_prefix}.bim"), path("${output_prefix}.fam")

    script:
        input_prefix = get_prefix(bed_file)
        output_prefix = "${input_prefix}.shared"
        """
        /opt/miniforge3/bin/mamba run -n plink2_env plink2 \
            --bfile ${input_prefix} \
            --extract range ${joint_variants} \
            --output-chr chr26 \
            --make-bed \
            --out ${output_prefix}
        """
}