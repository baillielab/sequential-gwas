process IndexReferenceGenome {
    publishDir "${params.GATK_PUBLISH_DIR}", mode: 'symlink'

    input:
        path reference_genome

    output:
        path("${reference_genome}.fai")

    script:
        """
        samtools faidx ${reference_genome}
        """
}