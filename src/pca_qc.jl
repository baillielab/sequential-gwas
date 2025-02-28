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

function get_outliers_and_save_plots(loadings; 
    output_prefix="before_pca_qc",
    iqr_factor=3,
    npcs=10
    )
    outliers = Set([])  
    fig = Figure(size=(1000, 1000))
    ax_row = 1
    ax_col = 1
    for pc in 1:npcs
        ax = Axis(fig[ax_row, ax_col],  ylabel="Loading", title="PC $pc")
        abs_loadings = abs.(loadings[!, "pc$(pc)_loading"])
        q1, q3 = quantile(abs_loadings, [0.25, 0.75])
        iqr_up_bound = (q3 + iqr_factor * (q3 - q1))
        loadings.IS_OUTLIER = abs_loadings .> iqr_up_bound
        union!(outliers, loadings.SNP[loadings.IS_OUTLIER])
        loadings.COLOR = ifelse.(loadings.IS_OUTLIER, "red", "blue")
        scatter!(ax, 1:nrow(loadings), abs_loadings, color=loadings.COLOR)
        ax_row, ax_col = isodd(pc) ? (ax_row, 2) : (ax_row+1, 1)
    end
    save(string(output_prefix, ".loadings.png"), fig)
    return outliers
end

plink2_pca(input_prefix, output_prefix; npcs=10) = run(Cmd([
    "plink2",
    "--bfile", input_prefix,
    "--pca", string(npcs),
    "--out", output_prefix])
)

function gcta_loadings(bfile_prefix, pca_prefix, outprefix)
    run(Cmd([
    "gcta64", 
    "--bfile", bfile_prefix, 
    "--pc-loading", pca_prefix, 
    "--out", outprefix
    ]))
    return CSV.read(string(outprefix, ".pcl"), DataFrame)
end

function plink2_exclude(outliers, input_prefix, outprefix)
    open("variants_to_exclude.txt", "w") do io
        for variant_id in outliers
            println(io, variant_id)
        end
    end
    run(Cmd([
        "plink2",
        "--bfile", input_prefix,
        "--exclude", "variants_to_exclude.txt",
        "--output-chr", "chr26",
        "--make-bed",
        "--out", outprefix
    ]))
end

function pca_qc(input_prefix, ancestry_file; 
    npcs=10, 
    iqr_factor=3,
    output_prefix = string(input_prefix, ".after_pca_qc")
    )
    # run PCA, compute PCA loadings and get outlier variants
    plink2_pca(input_prefix, input_prefix)
    loadings = gcta_loadings(input_prefix, input_prefix, input_prefix)
    outliers = get_outliers_and_save_plots(loadings; 
        output_prefix=string(input_prefix, ".before_pca_qc"), 
        iqr_factor=iqr_factor, 
        npcs=npcs
    )
    # Exclude outliers, run PCA and loadings and plot again to verify the effect of the operation
    plink2_exclude(outliers, input_prefix, output_prefix)
    plink2_pca(output_prefix, output_prefix)
    loadings = gcta_loadings(output_prefix, output_prefix, output_prefix)
    get_outliers_and_save_plots(loadings; 
        output_prefix=output_prefix,
        iqr_factor=iqr_factor, 
        npcs=npcs
    )
    # Further plots
    plot_pca(string(output_prefix, ".eigenvec"), ancestry_file; outprefix=output_prefix)
end