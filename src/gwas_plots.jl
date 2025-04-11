function harmonize(results)
    return DataFrames.select(results, 
        :CHROM => (x -> string.(x)) => :CHR,
        :GENPOS => :BP,
        :ID => :SNP,
        :LOG10P => (x -> exp10.(-x))  => :P,
    )
end

function gwas_plots(results_path, group; output_prefix = group)
    results = harmonize(CSV.read(results_path, DataFrame))
    # Plot Manhattan
    fig = Figure(size = (600, 400))
    ax = Axis(fig[1, 1], xlabel="Chromosome")
    GeneticsMakie.plotgwas!(ax, results, build=38)
    hidespines!(ax, :t, :r)
    Label(fig[1, 1, Top()], text = group, fontsize = 20)
    resize_to_layout!(fig)
    save(string(output_prefix, ".manhattan.png"), fig)
    # Plot QQ
    fig = Figure(size = (600, 600))
    ax = Axis(fig[1, 1])
    GeneticsMakie.plotqq!(ax, results; ystep = 5)
    hidespines!(ax, :t, :r)
    Label(fig[1, 1, Top()], text = group, fontsize = 20)
    resize_to_layout!(fig)
    save(string(output_prefix, ".qq.png"), fig)
end

function merge_regenie_chr_results(input_prefix; output="results.csv")
    dir, input_prefix_ = splitdir(input_prefix)
    dir = dir == "" ? "." : dir
    files = filter(x -> startswith(x, input_prefix_), readdir(dir))
    results = mapreduce(f -> CSV.read(joinpath(dir, f), DataFrame), vcat, files)
    CSV.write(output, results)
end