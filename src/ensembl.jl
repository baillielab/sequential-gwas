const ENSEMBL_SERVER = "https://rest.ensembl.org"

function get_vep(ids; ext="/vep/human/id")
    headers=Dict("Content-Type" => "application/json", "Accept" => "application/json")
    body = JSON.json(Dict("ids" => ids, "AlphaMissense" => 1, "Enformer" => 1))
    r = HTTP.post(ENSEMBL_SERVER*ext, headers, body)
    return JSON.parse(String(r.body))
end


function get_genomic_features(region; features=["gene", "transcript", "cds", "exon", "regulatory", "motif"])
    ext = string("/overlap/region/human/", region, "?", join(map(f -> "feature=$(f)", features), ";"))
    headers=Dict("Content-Type" => "application/json", "Accept" => "application/json")
    r = HTTP.get(ENSEMBL_SERVER*ext, headers)
    return JSON.parse(String(r.body))
end

function integrated_plot(variants_info, genomic_features)
    feature_colors = Dict("gene" => :blue, "transcript" => :green, "cds" => :orange, "exon" => :red, "regulatory" => :purple, "motif" => :brown)
    credible_sets = sort(unique(collect(variants_info.CS)))
    # Plot
    fig = Figure(size = (1000, 800))
    # GWAS P-values
    ax1 = Axis(fig[1, 1]; ylabel="-log10(P)", 
        xticksvisible=false, 
        xticklabelsvisible=false, 
        xgridvisible=false,
        ygridvisible=false
    )
    scatter!(ax1, 
        collect(variants_info.POS),
        collect(variants_info.LOG10P)
    )
    hlines!(ax1, 8, color=:green)
    # Fine Mapping PIPs
    credible_sets_colors = distinguishable_colors(length(credible_sets)+1, [RGB(1,1,1), RGB(0,0,0)])[2:end];
    cs_to_color = Dict(credible_sets .=> credible_sets_colors)
    ax2 = Axis(fig[2, 1]; ylabel="PIP", 
        xticksvisible=false, 
        xticklabelsvisible=false, 
        xgridvisible=false,
        ygridvisible=false
    )
    scatter!(ax2, 
        collect(variants_info.POS), 
        collect(variants_info.PIP); 
        color=[cs_to_color[cs] for cs in variants_info.CS], 
        markersize=10, 
        colormap=:tab10
    )
    pips_legend_elements = [(string("CS ", cs), PolyElement(color=cs_color, colorrange=1:10, colormap=:tab10)) for (cs, cs_color) in enumerate(credible_sets_colors[2:end])]
    Legend(fig[2, 2], getindex.(pips_legend_elements, 2), getindex.(pips_legend_elements, 1), "Fine Mapping CSs", framevisible = false)
    # Genomic annotations
    feature_types = unique(x["feature_type"] for x in genomic_features)
    feature_type_y_coord = Dict(feature_types .=> 1:length(feature_types))
    ax3 = Axis(fig[3, 1]; 
        xlabel="Position", 
        ylabel="Feature\nTrack",
        yticklabelsvisible=false,
        yticksvisible=false,
        xgridvisible=false,
        ygridvisible=false
    )
    for feature in genomic_features
        y_coord = feature_type_y_coord[feature["feature_type"]]
        lines!(ax3, 
            [feature["start"], feature["end"]], 
            [y_coord, y_coord], 
            color=feature_colors[feature["feature_type"]],
            linewidth=10
        )
    end
    ann_legend_elements = [(feature, LineElement(color=color)) for (feature, color) in feature_colors]
    Legend(fig[3, 2], getindex.(ann_legend_elements, 2), getindex.(ann_legend_elements, 1), "Genomic Features", framevisible = false)

    return fig
end

results_files = [
    "test/assets/EUR.SEVERE_COVID_19.rs7515509.variants_info.tsv",
    "test/assets/EUR.SEVERE_COVID_19.rs61738875.variants_info.tsv",
    "test/assets/EUR.SEVERE_COVID_19.rs114301457.variants_info.tsv",
    "test/assets/EUR.SEVERE_COVID_19.rs141942982.variants_info.tsv",
    "test/assets/EUR.SEVERE_COVID_19.rs189087412.variants_info.tsv"
]
features = ["gene", "exon", "regulatory", "motif"]

for result_file in results_files
    variants_info = CSV.read(result_file, DataFrame, delim=",")
    region = string(first(variants_info[!, "#CHROM"]), ":", minimum(variants_info.POS), "-", maximum(variants_info.POS))
    genomic_features = get_genomic_features(region; features=features)
    fig = integrated_plot(variants_info, genomic_features)
    display(fig)
end
