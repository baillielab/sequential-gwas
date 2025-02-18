function plot_pca(eigenvectors_file; outprefix="pca")
    eigenvectors = CSV.read(eigenvectors_file, DataFrame)
    # PC1 vs PC2 plot
    fig = Figure(size=(800, 800))
    ax = Axis(fig[1, 1], xlabel="PC1", ylabel="PC2")
    scatter!(ax, eigenvectors.PC1, eigenvectors.PC2)
    save(string(outprefix, ".1vs2.png"), fig)
    # All PCs plots
    n_pcs = size(eigenvectors, 2) - 2
    fig = Figure(size=(1000, 1000))
    ax_row = 1
    ax_col = 1
    for pc_k in 1:n_pcs-1
        ax = Axis(fig[ax_row, ax_col], xlabel="PC$(pc_k)", ylabel="PC$(pc_k+1)")
        scatter!(ax, eigenvectors[!, "PC$(pc_k)"], eigenvectors[!, "PC$(pc_k+1)"])
        ax_row, ax_col = if iseven(pc_k)
            (ax_row + 1, 1)
        else
            (ax_row, 2)
        end
    end
    save(string(outprefix, ".all.png"), fig)
end