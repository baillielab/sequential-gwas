function process_age(ages)
    return map(ages) do age
        if age == "NA"
            return missing
        else
            return parse(Int, age)
        end
    end
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

function read_and_process_covariates(covariates_file)
    covariates = CSV.read(
        covariates_file, 
        DataFrame
        )
    return select(covariates,
        :genotype_file_id => :IID,
        :age_years => process_age => :AGE,
        :sex => process_sexes => :SEX,
        :severe_cohort_primary_diagnosis => process_primary_diagnosis => :GENOMICC_PRIMARY_DIAGNOSIS,
        :cohort => process_cohort => :COHORT,
        :case_or_control => process_severity => :DIAGNOSIS_IS_SEVERE,
        :isaric_cohort_max_severity_score => process_isaric_score => :ISARIC_MAX_SEVERITY_SCORE
    )
end

function read_and_process_ancestry(ancestry_file)
    ancestries = CSV.read(ancestry_file, DataFrame)
    return select(ancestries, :FID, :IID, :Superpopulation => :ANCESTRY)
end

function read_and_process_pcs(pcs_file)
    return  CSV.read(pcs_file, DataFrame, drop=["#FID"])
end


function combine_covariates(
    covariates_file, 
    ancestry_file, 
    pcs_file; 
    output="covariates.merged.csv"
    )
    covariates = read_and_process_covariates(covariates_file)
    ancestries = read_and_process_ancestry(ancestry_file)
    pcs = read_and_process_pcs(pcs_file)
    merged = innerjoin(
        innerjoin(
            covariates, 
            ancestries, 
            on=:IID
        ),
        pcs,
        on=:IID
    )
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
    elseif row.COHORT == "uk"
        return 0
    else
        throw(ArgumentError("Unknown cohort"))
    end
end

function define_phenotype!(covariates, phenotype)
    if phenotype == "SEVERE_COVID_19"
        covariates.SEVERE_COVID_19 = map(eachrow(covariates)) do row
            is_severe_covid_19(row)
        end
    else
        throw(ArgumentError("Unsupported phenotype"))
    end
end

function process_covariates!(covariates, variables)
    covariate_processing_fns = Dict(
        "AGE" => v -> coalesce.(v, mean(skipmissing(v))),
        "SEX" => identity
    )
    if !issubset(variables, keys(covariate_processing_fns))
        throw(ArgumentError("Can only process the following covariates: ", keys(covariate_processing_fns)...))
    end
    for variable in variables
        covariates[!, variable] = covariate_processing_fns[variable](covariates[!, variable])
    end
end

function make_gwas_groups(covariates_file, variables_file; output_prefix="gwas", min_group_size=100)
    covariates = CSV.read(covariates_file, DataFrame)
    variables = YAML.load_file(variables_file)
    define_phenotype!(covariates, variables["phenotype"])
    process_covariates!(covariates, variables["covariates"])
    for (groupkey, group) in pairs(groupby(covariates, variables["group"], skipmissing=true))
        group_id = join(groupkey, "_")
        # Only retain individuals with no missing values for all covariates and phenotypes
        nomissing = dropmissing(select(group, "FID", "IID", variables["covariates"], variables["phenotype"]))
        if nrow(nomissing) < min_group_size
            @info "Skipping group $group_id because it has fewer than $min_group_size individuals."
            continue
        end
        # Write group covariates
        group_covariates = select(nomissing, "FID", "IID", variables["covariates"])
        CSV.write(string(output_prefix, ".covariates.", group_id, ".csv"), group_covariates)
        # Write group phenotype
        group_phenotype = select(nomissing, "FID", "IID", variables["phenotype"])
        CSV.write(string(output_prefix, ".phenotype.", group_id, ".csv"), group_phenotype)
        # Write group individuals
        CSV.write(string(output_prefix, ".individuals.", group_id, ".txt"), select(nomissing, ["FID", "IID"]), header=false, delim="\t")
    end
end