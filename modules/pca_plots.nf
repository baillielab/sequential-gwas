process PCAPlots {
    publishDir "${params.PCA_PUBLISH_DIR}/pca_plots", mode: 'symlink'

    input:
        tuple path(eigenvec), path(eigenval)
        path ancestry

    output:
        tuple path("pca.1vs2.png"), path("pca.all.png")

    script:
        """
        ${params.JULIA_CMD} plot-pca \
            ${eigenvec} \
            ${ancestry} \
            --outprefix pca
        """

}