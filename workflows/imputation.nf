include { get_prefix; get_julia_cmd } from '../modules/utils.nf'

process ImputeGenotypes {
    input:
        path topmed_api_token_file
        path genotypes

    script:
    jobname = get_prefix(genotypes[0])
    """
    TOPMED_API_TOKEN=`cat ${topmed_api_token_file}`
    curl https://imputation.biodatacatalyst.nhlbi.nih.gov/api/v2/jobs/submit/imputationserver \
        -X "POST" \
        -H "X-Auth-Token: \${TOPMED_API_TOKEN}" \
        -F "job-name=${jobname}" \
        -F "files=@${genotypes}" \
        -F "refpanel=apps@topmed-r3" \
        -F "build=hg38" \
        -F "phasing=eagle" \
        -F "population=all" \
        -F "meta=yes"
    """
}

process SplitByChromosome {
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
            --output-chr chr26 \
            --chr ${chr} \
            --recode vcf \
            --out ${output_prefix}
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

workflow Imputation {
    topmed_api_token = file("assets/topmed-api-token")
    bed_genotypes = Channel.fromPath("${params.GENOTYPES_PREFIX}.{bed,bim,fam}").collect()
    // chrs = GetChromosomes(bed_genotypes).splitText() { it.trim() }
    chrs = Channel.of(1..22).map {it -> "chr$it"}
    vcf_genotypes = SplitByChromosome(bed_genotypes, chrs)
    ImputeGenotypes(topmed_api_token, vcf_genotypes)
}