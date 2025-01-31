include { get_prefix} from './utils.nf'

process MakeSharedVariantsList {
    publishDir "results/array-genotypes/qced_shared_variants", mode: 'symlink'
    input:
        path genotypes_files

    output:
        path("variants_intersection.txt")

    script:
        """
        julia --project=/opt/sequential-gwas/ /opt/sequential-gwas/bin/seq-gwas.jl \
            write-variants-intersection
        """
}

process ExtractSharedVariants {
    publishDir "results/array-genotypes/qced_shared_variants", mode: 'symlink'

    input:
        tuple path(bed_file), path(bim_file), path(fam_file)
        path shared_variants

    output:
        tuple path("${output_prefix}.bed"), path("${output_prefix}.bim"), path("${output_prefix}.fam")

    script:
        input_prefix = get_prefix(bed_file)
        output_prefix = "${input_prefix}.shared"
        """
        /opt/miniforge3/bin/mamba run -n plink_env plink \
            --bfile ${input_prefix} \
            --noweb \
            --extract ${shared_variants} \
            --make-bed \
            --out ${output_prefix}
        """
}