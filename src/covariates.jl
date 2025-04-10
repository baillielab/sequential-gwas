function process_age(ages)
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

function process_sexes(sexes)
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

function read_and_process_covariates(covariates_file, inferred_covariates_file, required_covariates)
    # Read the covariates file
    covariates = CSV.read(
        covariates_file, 
        DataFrame
        )
    # Process the covariates
    select!(covariates,
        :genotype_file_id => :IID,
        :age_years => process_age => :AGE,
        :sex => process_sexes => :SEX,
        :severe_cohort_primary_diagnosis => process_primary_diagnosis => :GENOMICC_PRIMARY_DIAGNOSIS,
        :cohort => process_cohort => :COHORT,
        :case_or_control => process_severity => :DIAGNOSIS_IS_SEVERE,
        :isaric_cohort_max_severity_score => process_isaric_score => :ISARIC_MAX_SEVERITY_SCORE
    )
    # Read inferred covariates
    inferred_covariates = CSV.read(
        inferred_covariates_file, 
        DataFrame,
        select=[:FID, :IID, :ANCESTRY_ESTIMATE, :AFR, :SAS, :EAS, :AMR, :EUR, :PLATFORM]
    )
    # Join
    covariates = innerjoin(covariates, inferred_covariates, on=:IID)
    # Add user defined covariates
    add_user_defined_covariates!(covariates, required_covariates)

    return covariates
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
        for iid in SequentialGWAS.read_fam(file).IID
            sample_id_to_platform[iid] = "GSA-MD-24v3-0_A1"
        end
    end
    for iid in SequentialGWAS.read_fam(release_2024_now_fam).IID
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
    output="covariates.merged.csv"
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

function define_phenotypes!(covariates, phenotypes)
    for phenotype in phenotypes
        if phenotype == "SEVERE_COVID_19"
            covariates.SEVERE_COVID_19 = map(eachrow(covariates)) do row
                is_severe_covid_19(row)
            end
        else
            throw(ArgumentError("Unsupported phenotype"))
        end
    end
end

function add_user_defined_covariates!(covariates, variables)
    covariate_names = names(covariates)
    for variable in variables
        variable ∈ covariate_names && continue
        if occursin("_x_", variable)
            base_variables = split(variable, "_x_")
            issubset(base_variables, covariate_names) || throw(ArgumentError("Some base covariate in $variable was not found."))
            columns = (covariates[!, v] for v in base_variables)
            covariates[!, variable] = .*(columns...)
        else
            throw(ArgumentError("Covariate $variable not suported."))
        end
    end
end

function write_covariates_and_phenotypes_group(data, required_covariates, required_phenotypes; group_id="all", output_prefix="gwas", min_group_size=100)
    # Only retain individuals with no missing values for all covariates and phenotypes
    nomissing = dropmissing(select(data, "FID", "IID", required_covariates, required_phenotypes))
    if nrow(nomissing) < min_group_size
        @info "Skipping group $group_id because it has fewer than $min_group_size individuals."
        return
    end
    # Write group covariates
    covariates = select(nomissing, "FID", "IID", required_covariates)
    CSV.write(string(output_prefix, ".covariates.", group_id, ".csv"), covariates, delim="\t")
    # Write group phenotype
    phenotypes = select(nomissing, "FID", "IID", required_phenotypes)
    CSV.write(string(output_prefix, ".phenotype.", group_id, ".csv"), phenotypes, delim="\t")
    # Write group individuals
    CSV.write(string(output_prefix, ".individuals.", group_id, ".txt"), select(nomissing, ["FID", "IID"]), header=false, delim="\t")
end

function make_gwas_groups(
    covariates_file, 
    inferred_covariates_file, 
    variables_file; 
    output_prefix="gwas", 
    min_group_size=100
    )
    variables = YAML.load_file(variables_file)
    covariates = SequentialGWAS.read_and_process_covariates(covariates_file, inferred_covariates_file, variables["covariates"])
    SequentialGWAS.define_phenotypes!(covariates, variables["phenotypes"])
    if haskey(variables, "groupby")
        for (groupkey, group) in pairs(groupby(covariates, variables["groupby"], skipmissing=true, sort=true))
            group_id = join(groupkey, "_")
            write_covariates_and_phenotypes_group(group, variables["covariates"], variables["phenotypes"]; 
                group_id=group_id, 
                output_prefix=output_prefix, 
                min_group_size=min_group_size
            )
        end
    else
        write_covariates_and_phenotypes_group(covariates, variables["covariates"], variables["phenotypes"]; 
                group_id="all", 
                output_prefix=output_prefix, 
                min_group_size=min_group_size
        )
    end

end

function merge_covariates_pcs(covariates_file, pcs_file; output="covariates_pcs.csv")
    covariates = CSV.read(covariates_file, DataFrame)
    pcs = CSV.read(pcs_file, DataFrame, drop=["#FID"])
    merged = innerjoin(covariates, pcs, on=:IID)
    CSV.write(output, merged, delim="\t")
    return merged
end