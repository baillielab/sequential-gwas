function parse_pvalue(log10_pval::AbstractString)
    if log10_pval == "NA"
        return NaN
    else
        return parse_pvalue(parse(Float64, log10_pval))
    end
end

parse_pvalue(log10_pval::Real) = exp10(-log10_pval)

parse_a1freq(freq::AbstractString) = freq == "NA" ? NaN : parse(Float64, freq)

parse_a1freq(freq::Real) = freq

function harmonize(results; maf=0.01)
    harmonized_results =  DataFrames.select(results, 
        :CHROM => (x -> string.(x)) => :CHR,
        :GENPOS => :BP,
        :ID => :SNP,
        :LOG10P => (x -> parse_pvalue.(x))  => :P,
        :A1FREQ => (x -> parse_a1freq.(x)) => :A1FREQ
    )
    filtered_results = filter(
        x -> !isnan(x.P) && x.A1FREQ > maf, 
        harmonized_results
    )
    return filtered_results[!, [:CHR, :BP, :SNP, :P]]
end

function gwas_plots(results_path; maf=0.01, output_prefix = "gwas.plot")
    group_phenotype_string = replace(splitext(basename(results_path))[1], "regenie.results." => "")
    group, phenotype = split(group_phenotype_string, ".")
    results = harmonize(CSV.read(results_path, DataFrame), maf=maf)
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

function merge_regenie_chr_results(merge_list_file; output = "regenie.results.tsv")
    merge_list = readlines(merge_list_file)
    # Concatenating the files
    results = mapreduce(f -> CSV.read(f, DataFrame), vcat, merge_list)
    # Writing the output
    CSV.write(output, results; delim="\t", header=true)
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
    phenotypes_string="SEVERE_COVID_19",
    output_prefix="gwas", 
    min_cases_controls=100
    )
    phenotypes = split(phenotypes_string, ",")
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
    n_groups_passed = 0
    if groupby_string !== nothing
        groupby_variables = split(groupby_string, ",")
        for (groupkey, group) in pairs(groupby(covariates, groupby_variables, skipmissing=true, sort=true))
            group_id = join(groupkey, "_")
            n_phenotypes_passed = write_covariates_and_phenotypes_group(group, required_covariate_variables;
                group_id=group_id,
                phenotypes=phenotypes,
                output_prefix=output_prefix,
                min_cases_controls=min_cases_controls
            )
            n_groups_passed += n_phenotypes_passed
        end
    else
        n_groups_passed = write_covariates_and_phenotypes_group(covariates, required_covariate_variables; 
                group_id="all",
                phenotypes=phenotypes, 
                output_prefix=output_prefix, 
                min_cases_controls=min_cases_controls
        )
    end

    n_groups_passed > 0 || throw(ArgumentError("No group passed the min cases/controls threshold."))

    return 0
end
