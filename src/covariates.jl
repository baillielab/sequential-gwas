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
    return select(ancestries, :IID, :Superpopulation => :ANCESTRY)
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