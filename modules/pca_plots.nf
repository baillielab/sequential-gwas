process PCAPlots {
    publishDir "${params.PCA_PUBLISH_DIR}/pca_plots", mode: 'symlink'

    input:
        tuple path(eigenvec), path(eigenval)

    output:
        tuple path("pca.1vs2.png"), path("pca.all.png")

    script:
        """
        julia --project=/opt/sequential-gwas/ /opt/sequential-gwas/bin/seq-gwas.jl \
            plot-pca \
            ${eigenvec} \
            --outprefix pca
        """

}