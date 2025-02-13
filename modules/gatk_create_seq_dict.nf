process CreateSequenceDictionary {
    publishDir "${params.GATK_PUBLISH_DIR}", mode: 'symlink'

    input:
        path reference_genome

    output:
        path("${reference_genome.baseName}.dict")

    script:
        """
        gatk CreateSequenceDictionary -R ${reference_genome}
        """
}