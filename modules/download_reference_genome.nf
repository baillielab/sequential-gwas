process DownloadOrAccessReferenceGenome {
    storeDir "${params.GATK_DIR}"

    output:
        path "Homo_sapiens_assembly38.fasta"

    script:
        """
        wget -O Homo_sapiens_assembly38.fasta https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta
        """
}