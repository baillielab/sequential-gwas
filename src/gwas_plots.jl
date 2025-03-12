function harmonize(results)
    return select(results, 
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