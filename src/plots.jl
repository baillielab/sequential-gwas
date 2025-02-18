function plot_pc_k_vs_kplus1!(ax, eigenvectors, k, colormap)
    for (key, group) in pairs(groupby(eigenvectors, :Superpopulation))
        scatter!(
            ax, 
            group[!, "PC$(k)"],
            group[!, "PC$(k+1)"], 
            color=colormap[key.Superpopulation], 
            label=key.Superpopulation,
            colormap=:tab10,
            colorrange=(1, length(colormap))
        )
    end
end

function plot_pca(eigenvectors_file, ancestry_file; outprefix="pca")
    eigenvectors = CSV.read(eigenvectors_file, DataFrame)
    ancestry = CSV.read(ancestry_file, DataFrame, select=[:IID, :Superpopulation])
    leftjoin!(eigenvectors, ancestry, on=:IID)
    colormap = Dict(key => k for (k, key) in enumerate(unique(eigenvectors.Superpopulation)))
    # PC1 vs PC2 plot
    fig = Figure(size=(800, 800))
    ax = Axis(fig[1, 1], xlabel="PC1", ylabel="PC2")
    plot_pc_k_vs_kplus1!(ax, eigenvectors, 1, colormap)
    fig[2, :] = Legend(fig, ax, orientation=:horizontal, tellwidth=false)
    save(string(outprefix, ".1vs2.png"), fig)
    # All PCs plots
    n_pcs = size(eigenvectors, 2) - 3
    fig = Figure(size=(1000, 1000))
    ax_row = 1
    ax_col = 1
    local ax
    for pc_k in 1:n_pcs-1
        ax = Axis(fig[ax_row, ax_col], xlabel="PC$(pc_k)", ylabel="PC$(pc_k+1)")
        plot_pc_k_vs_kplus1!(ax, eigenvectors, pc_k, colormap)
        ax_row, ax_col = if iseven(pc_k)
            (ax_row + 1, 1)
        else
            (ax_row, 2)
        end
    end
    fig[ax_row+1, :] = Legend(fig, ax, 
        orientation=:horizontal,
        tellwidth = false)
    save(string(outprefix, ".all.png"), fig)
end