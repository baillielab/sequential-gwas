function merge_chr_results(gwas_merge_list_file, finemapping_merge_list_file; output_prefix = "results.all_chr")
    # Concatenate GWAS files
    gwas_merge_list = readlines(gwas_merge_list_file)
    results = mapreduce(f -> CSV.read(f, DataFrame), vcat, gwas_merge_list)
    CSV.write(string(output_prefix, ".gwas.tsv"), results; delim="\t", header=true)
    # Concatenate finemapping files
    finemapping_merge_list = readlines(finemapping_merge_list_file)
    finemapping_results = mapreduce(f -> CSV.read(f, DataFrame), vcat, finemapping_merge_list)
    CSV.write(string(output_prefix, ".finemapping.tsv"), finemapping_results; delim="\t", header=true)

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
