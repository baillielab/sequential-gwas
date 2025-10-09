process WGSIndividuals  {
    publishDir "${params.WGS_PUBLISH_DIR}", mode: 'symlink'

    input:
        path(gvcf_files)

    output:
        path("gvcf_sample_ids.txt")

    script:
        """
        for file in *.gvcf.gz; do
            bcftools query -l "\$file"
        done | sort | uniq > gvcf_sample_ids.txt
        """
}