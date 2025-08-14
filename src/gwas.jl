function parse_pvalue(log10_pval::AbstractString)
    if log10_pval == "NA"
        return NaN
    else
        return parse_pvalue(parse(Float64, log10_pval))
    end
end

parse_pvalue(log10_pval::Real) = exp10(-log10_pval)

function harmonize(results)
    harmonized_results =  DataFrames.select(results, 
        :CHROM => (x -> string.(x)) => :CHR,
        :GENPOS => :BP,
        :ID => :SNP,
        :LOG10P => (x -> parse_pvalue.(x))  => :P,
    )
    return filter(:P => !(isnan), harmonized_results)
end

function gwas_plots(results_path; output_prefix = "gwas.plot")
    group_phenotype_string = replace(splitext(basename(results_path))[1], "regenie.results." => "")
    group, phenotype = split(group_phenotype_string, ".")
    results = harmonize(CSV.read(results_path, DataFrame))
    # Plot Manhattan
    fig = Figure(size = (600, 400))
    ax = Axis(fig[1, 1], xlabel="Chromosome")
    GeneticsMakie.plotgwas!(ax, results, build=38)
    hidespines!(ax, :t, :r)
    Label(fig[1, 1, Top()], text = group, fontsize = 20)
    resize_to_layout!(fig)
    save(string(output_prefix, ".$group.$phenotype.manhattan.png"), fig)
    # Plot QQ
    fig = Figure(size = (600, 600))
    ax = Axis(fig[1, 1])
    GeneticsMakie.plotqq!(ax, results; ystep = 5)
    hidespines!(ax, :t, :r)
    Label(fig[1, 1, Top()], text = group, fontsize = 20)
    resize_to_layout!(fig)
    save(string(output_prefix, ".$group.$phenotype.qq.png"), fig)

    return 0
end

function merge_regenie_chr_results(merge_list_file; output_prefix = "regenie.results")
    merge_list = CSV.read(merge_list_file, DataFrame; header=["FILE"])
    merge_list.BASENAME = basename.(merge_list.FILE)
    merge_list.GROUP = getindex.(split.(merge_list.BASENAME, "."), 1)
    merge_list.PHENOTYPE = replace.(getindex.(split.(merge_list.BASENAME, "."), 3), "step2_" => "")
    for (groupkey, group) in pairs(groupby(merge_list, [:GROUP, :PHENOTYPE]))
        results = mapreduce(f -> CSV.read(f, DataFrame), vcat, group.FILE)
        output_file = string(output_prefix, ".", groupkey.GROUP, ".", groupkey.PHENOTYPE, ".tsv")
        CSV.write(output_file, results; delim="\t", header=true)
    end
    return 0
end

get_chr_out_string(pc_filename) = splitext(splitext(pc_filename)[1])[2][2:end]

function read_loco_pcs(pc_file)
    chr_out = GenomiccWorkflows.get_chr_out_string(pc_file)
    pcs = CSV.read(pc_file, DataFrame, drop=["#FID"])
    PC_colnames = filter(!=("IID"), names(pcs))
    for PC_colname in PC_colnames
        rename!(pcs, Symbol(PC_colname) => Symbol(string(uppercase(chr_out), "_", PC_colname)))
    end
    return pcs
end

function merge_covariates_and_pcs(covariates_file, pcs_prefix; output="covariates_and_pcs.csv")
    covariates = CSV.read(covariates_file, DataFrame)
    pcs_dir = dirname(pcs_prefix)
    pcs_dir, pcs_prefix = pcs_dir == "" ? (".", "./$pcs_prefix") : (pcs_dir, pcs_prefix)
    pcs_files = filter(startswith(pcs_prefix), readdir(pcs_dir, join=true))
    chrs = unique(GenomiccWorkflows.get_chr_out_string.(pcs_files))
    chr_out_pcs = map(chrs) do chr
        chr_out_pcs_files = filter(x -> occursin(chr, x), pcs_files)
        mapreduce(read_loco_pcs, vcat, chr_out_pcs_files)
    end
    covariates_and_pcs = innerjoin(covariates, chr_out_pcs..., on=:IID)
    CSV.write(output, covariates_and_pcs, delim="\t", missingstring="NA")
    return 0
end

function make_gwas_groups(
    covariates_file; 
    groupby_string=nothing,
    covariates_string="AGE",
    output_prefix="gwas", 
    min_group_size=100
    )

    # Define additional covariates
    covariates, required_covariate_variables = read_and_process_covariates(covariates_file; covariates_string=covariates_string)
    # Write new covariates to file
    CSV.write(
        string(output_prefix, ".covariates.csv"), 
        covariates, 
        delim="\t"
    )
    # Write required covariates list to file for REGENIE
    open(string(output_prefix, ".covariates_list.txt"), "w") do io
        for covariate in required_covariate_variables
            println(io, covariate)
        end
    end
    # Make groups
    if groupby_string !== nothing
        groupby_variables = split(groupby_string, ",")
        for (groupkey, group) in pairs(groupby(covariates, groupby_variables, skipmissing=true, sort=true))
            group_id = join(groupkey, "_")
            write_covariates_and_phenotypes_group(group; 
                group_id=group_id, 
                output_prefix=output_prefix, 
                min_group_size=min_group_size
            )
        end
    else
        write_covariates_and_phenotypes_group(covariates; 
                group_id="all", 
                output_prefix=output_prefix, 
                min_group_size=min_group_size
        )
    end
end
