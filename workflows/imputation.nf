include { get_prefix; get_julia_cmd } from '../modules/utils.nf'

process ImputeGenotypes {
    input:
        path topmed_api_token_file
        path genotypes

    script:
    """
    topmed_api_token=`cat ${topmed_api_token_file}`
    curl https://imputation.biodatacatalyst.nhlbi.nih.gov/api/v2/jobs/submit/imputationserver \
        -X "POST" \
        -H "X-Auth-Token: \${topmed_api_token}" \
        -F "job-name=Attempt 1" \
        -F "files=@${genotypes}" \
        -F "refpanel=apps@topmed-r3" \
        -F "build=hg38" \
        -F "phasing=eagle" \
        -F "population=all" \
        -F "meta=yes"
    """
}

process BedToChrVCF {
    input:
        path genotypes
        val chr

    output:
        path "${output_prefix}.vcf.gz"

    script:
        input_prefix = get_prefix(genotypes[0])
        output_prefix = "${input_prefix}.${chr}"
        """
        plink \
            --bfile ${input_prefix} \
            --chr ${chr} \
            --recode vcf \
            --out ${output_prefix}
        bgzip ${output_prefix}.vcf
        """
}

process GetChromosomes {
    input:
        path genotypes

    output:
        path "chromosomes.txt"

    script:
        input_prefix = get_prefix(genotypes[0])
        """
        ${get_julia_cmd(task.cpus)} write-chromosomes ${input_prefix}
        """
}

workflow ImputationWorkflow {
    topmed_api_token = file("assets/topmed-api-token")
    bed_genotypes = Channel.fromPath("results/array-genotypes/merged_qced/genotypes.merged.qced*").collect()
    chrs = GetChromosomes(bed_genotypes).splitText() { it.trim() }
    vcf_genotypes = BedToChrVCF(bed_genotypes, chrs)
    ImputeGenotypes(topmed_api_token, vcf_genotypes)
}