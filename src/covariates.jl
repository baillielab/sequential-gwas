function process_genomicc_sexes(sexes)
    return map(sexes) do sex
        if sex == "Male"
            return 1
        elseif sex == "Female"
            return 0
        else
            return missing
        end
    end
end

function process_ukb_age(years_of_birth)
    current_year = year(Dates.now())
    ages = current_year .- years_of_birth
    μ = round(Int, mean(skipmissing(ages)))
    return coalesce.(ages, μ)
end

function read_and_process_covariates(covariates_file;
    covariates_string=nothing
    )
    # Read the covariates file
    covariates = CSV.read(covariates_file, DataFrame)
    # Add user defined covariates
    required_covariate_variables = add_user_defined_covariates!(covariates, covariates_string)

    return covariates, required_covariate_variables
end

function read_and_process_ancestry(ancestry_file)
    ancestries = CSV.read(ancestry_file, DataFrame)
    rename!(ancestries, :Superpopulation => :ANCESTRY_ESTIMATE)
    return ancestries
end

function read_and_process_pcs(pcs_file)
    return CSV.read(pcs_file, DataFrame, drop=["#FID"])
end

function map_sample_ids_to_platform(wgs_samples_file,
    release_r8_fam,
    release_2021_2023_fam,
    release_2024_now_fam)
    # The order in which the files are processed is important
    # because platforms are given the following priority: WGS > GSA-MD-48v4-0_A1 > GSA-MD-24v3-0_A1 
    sample_id_to_platform = Dict()
    for file in (release_r8_fam, release_2021_2023_fam)
        for iid in GenomiccWorkflows.read_fam(file).IID
            sample_id_to_platform[iid] = "GSA-MD-24v3-0_A1"
        end
    end
    for iid in GenomiccWorkflows.read_fam(release_2024_now_fam).IID
        sample_id_to_platform[iid] = "GSA-MD-48v4-0_A1"
    end
    for iid in readlines(wgs_samples_file)
        sample_id_to_platform[iid] = "WGS"
    end
    return sample_id_to_platform
end

function combine_covariates(
    ancestry_file, 
    pcs_file,
    wgs_samples_file,
    release_r8_fam,
    release_2021_2023_fam,
    release_2024_now_fam;
    output="covariates.inferred.csv"
    )
    sample_id_to_platform = map_sample_ids_to_platform(
        wgs_samples_file,
        release_r8_fam,
        release_2021_2023_fam,
        release_2024_now_fam
    )
    ancestries = read_and_process_ancestry(ancestry_file)
    pcs = read_and_process_pcs(pcs_file)
    merged = innerjoin(
        ancestries,
        pcs,
        on=:IID
    )
    merged.PLATFORM = map(merged.IID) do iid
        sample_id_to_platform[iid]
    end
    CSV.write(output, merged)
    return merged
end

function add_user_defined_covariates!(covariates, covariates_string)
    required_covariate_variables = split(covariates_string, ",")
    all_colnames = names(covariates)
    updated_required_covariate_variables = []
    for variable in required_covariate_variables
        if variable ∈ all_colnames
            if eltype(covariates[!, variable]) <: AbstractString
                covariates[!, variable] = categorical(covariates[!, variable])
                mach = machine(OneHotEncoder(), covariates[!, [variable]])
                fit!(mach, verbosity=0)
                Xt = MLJBase.transform(mach)
                for colname in names(Xt)
                    covariates[!, colname] = Xt[!, colname]
                    push!(updated_required_covariate_variables, colname)
                end
            else
                push!(updated_required_covariate_variables, variable)
            end
        elseif occursin("_x_", variable)
            base_variables = split(variable, "_x_")
            issubset(base_variables, all_colnames) || throw(ArgumentError("Some base covariate in $variable was not found."))
            columns = (covariates[!, v] for v in base_variables)
            covariates[!, variable] = .*(columns...)
            push!(updated_required_covariate_variables, variable)
        else
            throw(ArgumentError("Covariate $variable not suported."))
        end
    end
    return updated_required_covariate_variables
end

function write_covariates_and_phenotypes_group(data; group_id="all", phenotypes=["SEVERE_COVID_19"], output_prefix="gwas", min_cases_controls=100)
    group_phenotypes = []
    for phenotype in phenotypes
        group_cases_controls_df = combine(groupby(data, phenotype, skipmissing=true), nrow)
        group_cases_controls_dict = Dict(
            string(val) => n for (val, n) in 
                zip(group_cases_controls_df[!, phenotype], group_cases_controls_df.nrow)
        ) # Stringify the phenotype values, in case the input is already in format [0, 1, NA]
        cases = get(group_cases_controls_dict, "1", 0)
        controls = get(group_cases_controls_dict, "0", 0)
        if controls < min_cases_controls || cases < min_cases_controls
            @info "Skipping phenotype $phenotype for group $group_id because it has fewer than $min_cases_controls cases/controls: (cases: $(cases), controls: $(controls))."
        else
            push!(group_phenotypes, phenotype)
        end
    end
    # Only retain individuals with no missing values for all covariates and phenotypes
    if isempty(group_phenotypes)
        @info "No phenotype passing the min cases/controls threshold for group $group_id, which will be skipped ."
        return
    end
    # Write group individuals
    CSV.write(string(output_prefix, ".individuals.", group_id, ".txt"), DataFrames.select(data, ["FID", "IID"]), header=false, delim="\t")
    # Write group phenotypes
    open(string(output_prefix, ".phenotypes.", group_id, ".txt"), "w") do io
        for phenotype in group_phenotypes
            println(io, phenotype)
        end
    end
end
