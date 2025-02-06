process CreateSequenceDictionary {
    publishDir "results/reference-genome", mode: 'symlink'

    input:
        path reference_genome

    output:
        path("${reference_genome.baseName}.dict")

    script:
        """
        gatk CreateSequenceDictionary -R ${reference_genome}
        """
}