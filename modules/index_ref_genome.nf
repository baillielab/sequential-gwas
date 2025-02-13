process IndexReferenceGenome {
    publishDir "${params.GATK_PUBLISH_DIR}", mode: 'symlink'

    input:
        path reference_genome

    output:
        path("${reference_genome}.fai")

    script:
        """
        /opt/miniforge3/bin/mamba run -n bcftools_env samtools faidx ${reference_genome}
        """
}