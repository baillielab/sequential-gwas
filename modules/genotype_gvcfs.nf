include { get_julia_cmd } from './utils.nf'

process GenotypeGVCFs{
    label "multithreaded"
    publishDir "${params.WGS_PUBLISH_DIR}/genotyped", mode: 'symlink'

    input:
        tuple val(prefix), path(gvcf)
        path joint_variants_gatk
        path shared_variants_plink
        path reference_genome
        path reference_genome_index
        path reference_genome_dict

    output:
        tuple path("${output_prefix}.bed"), path("${output_prefix}.bim"), path("${output_prefix}.fam")

    script:
        output_prefix = "${prefix}.shared"
        """
        # Genotype variants using GATK
        gatk GenotypeGVCFs \
            -R ${reference_genome} \
            -V ${prefix}.gvcf.gz \
            --intervals ${joint_variants_gatk} \
            --interval-padding 10 \
            --force-output-intervals ${joint_variants_gatk} \
            -O ${output_prefix}.vcf.gz
        
        # Converts to plink format
        plink2 \
            --threads ${task.cpus} \
            --memory ${task.memory.toMega().toString()} \
            --vcf ${prefix}.shared.vcf.gz \
            --output-chr chr26 \
            --max-alleles 2 \
            --set-all-var-ids @:# \
            --make-bed \
            --out ${output_prefix}.temp
        
        # Fills in missing allele and update variant ids
        ${get_julia_cmd(task.cpus)} complete-bim-with-ref \
            ${output_prefix}.temp.bim \
            ${shared_variants_plink} \
            --output ${output_prefix}.temp.bim

        # Only keep variants that are in the shared list
        plink2 \
            --threads ${task.cpus} \
            --memory ${task.memory.toMega().toString()} \
            --bfile ${output_prefix}.temp \
            --max-alleles 2 \
            --output-chr chr26 \
            --extract ${shared_variants_plink} \
            --make-bed \
            --out ${output_prefix}

        # Cleanup
        rm ${output_prefix}.temp.*
        rm ${output_prefix}.vcf.*
        """
}