function process_genomicc_age(ages)
    int_ages = map(ages) do age
        if age == "NA"
            return missing
        else
            return parse(Int, age)
        end
    end
    μ = round(Int, mean(skipmissing(int_ages)))
    return coalesce.(int_ages, μ)
end

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

function process_severity(severities)
    return map(severities) do severity
        if severity == "case"
            return 1
        elseif severity == "control"
            return 0
        else
            return missing
        end
    end
end

function process_primary_diagnosis(primary_diagnoses)
    return map(primary_diagnoses) do primary_diagnosis
        if primary_diagnosis == "NA"
            return missing
        else
            return primary_diagnosis
        end
    end
end

function process_isaric_score(isaric_scores)
    return map(isaric_scores) do isaric_score
        if isaric_score == "NA"
            return missing
        else
            return parse(Int, isaric_score)
        end
    end
end

function process_cohort(cohorts)
    replace(cohorts, "NA" => missing, "severe" => "genomicc")
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
    return  CSV.read(pcs_file, DataFrame, drop=["#FID"])
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

function is_severe_covid_19(row)
    if row.COHORT == "genomicc"
        if row.GENOMICC_PRIMARY_DIAGNOSIS == "COVID-19"
            return 1
        else
            return missing
        end
    elseif row.COHORT == "isaric4c"
        if row.ISARIC_MAX_SEVERITY_SCORE >= 4
            return 1
        else
            return missing
        end
    elseif row.COHORT == "mild"
        return 0
    elseif row.COHORT == "react"
        return 0
    elseif row.COHORT == "ukb"
        return 0
    else
        throw(ArgumentError("Unknown cohort"))
    end
end

function is_case(x)
    if x == "case"
        return 1
    elseif x == "control"
        return 0
    else
        return missing
    end
end

function define_phenotypes!(covariates, phenotypes)
    for phenotype in phenotypes
        if phenotype == "SEVERE_COVID_19"
            covariates.SEVERE_COVID_19 = map(eachrow(covariates)) do row
                is_severe_covid_19(row)
            end
        elseif phenotype == "case_or_control"
            covariates.case_or_control = map(covariates.case_or_control) do x
                is_case(x)
            end
        else
            throw(ArgumentError("Unsupported phenotype: $phenotype"))
        end
    end
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

function write_covariates_and_phenotypes_group(data; group_id="all", output_prefix="gwas", min_group_size=100)
    # Only retain individuals with no missing values for all covariates and phenotypes
    if nrow(data) < min_group_size
        @info "Skipping group $group_id because it has fewer than $min_group_size individuals."
        return
    end
    # Write group individuals
    CSV.write(string(output_prefix, ".individuals.", group_id, ".txt"), DataFrames.select(data, ["FID", "IID"]), header=false, delim="\t")
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


function read_loco_pcs(pc_file)
    chr_out = first(split(pc_file, "."))
    pcs = CSV.read(pc_file, DataFrame, drop=["#FID"])
    PC_colnames = filter(!=("IID"), names(pcs))
    for PC_colname in PC_colnames
        rename!(pcs, Symbol(PC_colname) => Symbol(string(uppercase(chr_out), "_OUT_", PC_colname)))
    end
    return pcs
end

function merge_covariates_and_pcs(covariates_file, pcs_prefix; output="covariates_and_pcs.csv")
    covariates = CSV.read(covariates_file, DataFrame)
    pcs_dir = dirname(pcs_prefix)
    pcs_dir = pcs_dir == "" ? "." : pcs_dir
    pcs_files = filter(startswith(pcs_prefix), readdir(pcs_dir))
    chrs = unique(getindex.(split.(pcs_files, "."), 1))
    chr_out_pcs = map(chrs) do chr
        chr_out_pcs_files = filter(x -> startswith(x, chr*"."), pcs_files)
        mapreduce(read_loco_pcs, vcat, chr_out_pcs_files)
    end
    covariates_and_pcs = innerjoin(covariates, chr_out_pcs..., on=:IID)
    CSV.write(output, covariates_and_pcs, delim="\t")
    return 0
end