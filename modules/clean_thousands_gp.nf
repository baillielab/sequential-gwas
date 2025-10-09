include { get_prefix } from './utils.nf'

process CleanKGP {
    label "multithreaded"
    
    input:
        tuple path(vcf_file), path(vcf_index)

    output:
        tuple path("${input_prefix}.clean.vcf.gz"), path("${input_prefix}.clean.vcf.gz.tbi")

    script:
        // (i) Splits multi-allelic variants into multiple rows 
        // (ii) resets the ID column to the CHROM:POS:REF:ALT format 
        // (iii) normalizes the VCF file
        // (iv) indexes the output
        input_prefix = get_prefix(vcf_file)[0..-5]
        output = "${input_prefix}.clean.vcf.gz"
        """
        bcftools norm --threads ${task.cpus} -m +any ${vcf_file} | \
        bcftools annotate --threads ${task.cpus} -x ID -I +'%CHROM:%POS:%REF:%ALT' | \
        bcftools norm --threads ${task.cpus} -Oz --rm-dup both > ${output}
        
        bcftools index --threads ${task.cpus} --tbi ${output}
        """
}