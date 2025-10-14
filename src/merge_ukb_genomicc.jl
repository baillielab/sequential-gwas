function get_severe_infections_map()
    severe_infections_map = Dict()
    for (raw_name, clean_name) in [
        "COVID_19" => "SEVERE_COVID_19",
        "INFLUENZA" => "SEVERE_INFLUENZA",
        "PNEUMONIA" => "SEVERE_PNEUMONIA",
        "PANCREATITIS" => "SEVERE_PANCREATITIS",
        "RSV" => "SEVERE_RSV",
        "SOFT_TISSUE_INFECTION" => "SEVERE_SOFT_TISSUE_INFECTION",
        "ECLS" => "SEVERE_ECLS",
        "REACTION_TO_VACCINATION" => "SEVERE_REACTION_TO_VACCINATION"
        ]
        severe_infections_map[raw_name] = (row -> is_severe_infection(row, raw_name), clean_name)
    end
    return severe_infections_map
end

format_chromosome!(bim) =
    bim.CHR_CODE = replace.(string.(bim.CHR_CODE), "chr" => "")

"""
    update_variant_ids_with_map!(bim, variant_ids_map)

Updates the variant IDs to match the KGP if the following conditions are met:
    - The position is present in the KGP dataset.
    - The UKB alleles match the KGP alleles.

These variant IDs are marked for deletion and returned as a set.
"""
function update_variant_ids_with_map!(bim, variant_ids_map)
    unmapped_ids = Set()
    bim.VARIANT_ID = map(zip(bim.CHR_CODE, bim.BP_COORD, bim.VARIANT_ID, bim.ALLELE_1, bim.ALLELE_2)) do (chr, loc, old_id, all_1, all_2)
        # If the variant's position is in the KGP
        if haskey(variant_ids_map, (chr, loc))
            new_variant_id, kgp_alleles = variant_ids_map[(chr, loc)]
            # If the alleles match, update the variant ID
            if Set([all_1, all_2]) == kgp_alleles
                new_variant_id
            # Otherwise keep the old ID but mark it for deletion
            else
                push!(unmapped_ids, old_id)
                old_id
            end
        # In other cases, keep the old ID and mark it for deletion
        else
            push!(unmapped_ids, old_id)
            old_id
        end
    end
    return unmapped_ids
end

function align_ukb_variants_with_kgp_and_keep_unrelated(ukb_bed_prefix, kgp_bed_prefix; out_prefix="ukb_unrelated", relatedness_degree=3)
    tmpdir = mktempdir()
    # Load KGP variants info
    kgp_bim = GenomiccWorkflows.read_bim(string(kgp_bed_prefix, ".bim"))
    format_chromosome!(kgp_bim)

    # Create a map of variant IDs from KGP
    kgp_variant_ids_map = Dict((chr, loc) => (v_id, Set([all_1, all_2])) for (chr, loc, v_id, all_1, all_2) in 
        zip(kgp_bim.CHR_CODE, kgp_bim.BP_COORD, kgp_bim.VARIANT_ID, kgp_bim.ALLELE_1, kgp_bim.ALLELE_2))

    # Load UKB variants info
    ukb_bim = GenomiccWorkflows.read_bim(string(ukb_bed_prefix, ".bim"))
    format_chromosome!(ukb_bim)

    # Map variant IDs to KGP if possible, otherwise they will be dropped
    unmapped_ids = update_variant_ids_with_map!(ukb_bim, kgp_variant_ids_map)

    # Find multi-allelic variants (split on multiple lines)
    multi_allelic_variants_df = filter(
        :nrow => >(1), 
        combine(groupby(ukb_bim, [:CHR_CODE, :BP_COORD]), nrow, :VARIANT_ID)
    )
    multi_allelic_variants = unique(multi_allelic_variants_df.VARIANT_ID)

    # Write variants to drop
    variants_to_drop_file = joinpath(tmpdir, "variants_to_drop.txt")
    variants_to_drop = collect(union(unmapped_ids, multi_allelic_variants))
    CSV.write(
        variants_to_drop_file, 
        DataFrame(VARIANT_ID=variants_to_drop),
        header=false
    )
    # Write new bim file
    new_bim_file = joinpath(tmpdir, "new.bim")
    CSV.write(
        new_bim_file,
        ukb_bim, 
        header=false, 
        delim="\t"
    )
    # Drop related individual using king
    run(`king --cpus $(nthreads()) -b $ukb_bed_prefix.bed --bim $new_bim_file --fam $ukb_bed_prefix.fam --unrelated --degree $relatedness_degree`)
    # Drop variants using plink2
    run(`plink2 --bed $ukb_bed_prefix.bed \
    --bim $new_bim_file \
    --fam $ukb_bed_prefix.fam \
    --exclude $variants_to_drop_file \
    --keep kingunrelated.txt \
    --output-chr chr26 \
    --make-bed \
    --out $out_prefix`
    )
    return 0
end


"""
    is_severe_infection(row, infection_name)

An individual is considered severe for an infection if it had the infection and is severely ill. 
If an individual does not have the infection, it is filled missing not to contaminate downstream analyses. 
They could be severe for another infection for example and are not appropriate controls.
"""
function is_severe_infection(row, infection_name)
    if row.PRIMARY_DIAGNOSIS == infection_name
        return row.IS_SEVERELY_ILL
    else
        return missing
    end
end

function is_severely_ill(cohort, severity_score)
    if cohort in ("GENOMICC_SEVERE", "GEN_INT_PAKISTAN")
        return 1
    elseif cohort == "ISARIC4C"
        if severity_score == "NA"
            return missing
        else
            score = parse(Int, severity_score)
            return score >= 4 ? 1 : 0
        end
    elseif cohort in ("GENOMICC_MILD", "GENOMICC_REACT")
        return 0
    else
        throw(ArgumentError("Unknown cohort: $cohort"))
    end
end

function add_is_severely_ill_col!(df)
    return transform!(df, 
        [:COHORT, :ISARIC_MAX_SEVERITY_SCORE] => ByRow(is_severely_ill) => :IS_SEVERELY_ILL
    )
end

"""
    alive_at_assessment(alive_60_genomicc, alive_28_isaric)

The provided column names are not ideal:
    - ALIVE_AT_28_DAYS: only concerns individuals with COHORT = ISARIC4C
    - ALIVE_AT_60_DAYS: only concerns individuals with COHORT = GenOMICC

The retention time is not the same for the two cohorts. 
The created column ALIVE_AT_ASSESSMENT informs wether an individual had died at the end of their respective assessment time (28 or 60 days).
"""
function alive_at_assessment(alive_60_genomicc, alive_28_isaric)
    alive_60_genomicc = lowercase(alive_60_genomicc)
    alive_28_isaric = lowercase(alive_28_isaric)

    if alive_60_genomicc == "yes"
        return 1
    elseif alive_60_genomicc == "no"
        return 0
    elseif alive_60_genomicc == "na"
        if alive_28_isaric == "yes"
            return 1
        elseif alive_28_isaric == "no"
            return 0
        elseif alive_28_isaric == "na"
            return missing
        else
            throw(ArgumentError("Unknown value for ALIVE_AT_28_DAYS: $alive_28_isaric"))
        end
    else
        throw(ArgumentError("Unknown value for ALIVE_AT_60_DAYS: $alive_60_genomicc"))
    end
end

function add_alive_at_assessment_col!(df)
    return transform!(df, 
        [:ALIVE_AT_60_DAYS, :ALIVE_AT_28_DAYS] => ByRow(alive_at_assessment) => :ALIVE_AT_ASSESSMENT
    )
end

function concat_ukb_covariates(ukb_covariates_file, ukb_inferred_covariates_file, genomicc_covariates, infection_names)
    ukb_covariates = CSV.read(ukb_covariates_file, DataFrame)
    ukb_inferred_covariates = CSV.read(ukb_inferred_covariates_file, DataFrame)
    ukb_all_covariates = innerjoin(
        ukb_covariates,
        ukb_inferred_covariates,
        on=:eid => :IID
    )
    n_ukb_samples = nrow(ukb_all_covariates)
    # Select aand process AGE, SEX and inferred ancestry
    DataFrames.select!(ukb_all_covariates,
        :eid => :FID,
        :eid => :IID,
        Symbol("34-0.0") => process_ukb_age => :AGE,
        Symbol("22001-0.0") => :SEX,
        :Superpopulation => :SUPERPOPULATION,
        :AFR,
        :AMR,
        :EAS,
        :EUR,
        :SAS,
    )
    # Severe infections are all 0 for UKB
    for infection_name in infection_names
        ukb_all_covariates[!, infection_name] = zeros(Int, n_ukb_samples)
    end
    # Fill cohort
    ukb_all_covariates.COHORT = fill("UKBIOBANK", n_ukb_samples)
    # Add ALIVE_AT_ASSESSMENT: we fill with missings
    if hasproperty(genomicc_covariates, :ALIVE_AT_ASSESSMENT)
        ukb_all_covariates.ALIVE_AT_ASSESSMENT = fill(missing, n_ukb_samples)
    end
    # Add PRIMARY_DIAGNOSIS
    ukb_all_covariates.PRIMARY_DIAGNOSIS = fill(missing, n_ukb_samples)
    # ADD IS_SEVERELY_ILL
    ukb_all_covariates.IS_SEVERELY_ILL = fill(0, n_ukb_samples)

    return vcat(genomicc_covariates, ukb_all_covariates)
end

function add_primary_diagnosis!(df)
    prim_diag_map = Dict(
        "covid-19" => "COVID_19",
        "isaric4c covid-19" => "COVID_19",
        "mild covid-19" => "COVID_19",
        "react covid-19" => "COVID_19",
        "influenza virus" => "INFLUENZA",
        "pneumonia with radiographic changes at presentation to critical care" => "PNEUMONIA",
        "pancreatitis of any aetiology" => "PANCREATITIS",
        "rsv (respiratory syncytial virus) infection" => "RSV",
        "soft tissue infections causing systemic sepsis" => "SOFT_TISSUE_INFECTION",
        "ecls" => "ECLS",
        "reaction to vaccination" => "REACTION_TO_VACCINATION"
    )
    df.PRIMARY_DIAGNOSIS = map(d -> prim_diag_map[d], df.PRIM_DIAGNOSIS_ODAP)
    return df
end

function clean_cohort!(df)
    df.COHORT = uppercase.(df.COHORT)
    return df
end

function process_genomicc_covariates(
    genomicc_covariates_file;
    ukb_covariates_file=nothing,
    ukb_inferred_covariates_file=nothing,
    output_file="covariates.processed.csv"
    )
    severe_infections_map = get_severe_infections_map()
    infection_names = last.(values(severe_infections_map))
    # Process GenOMICC covariates
    genomicc_covariates = CSV.read(genomicc_covariates_file, DataFrame)
    clean_cohort!(genomicc_covariates)
    add_is_severely_ill_col!(genomicc_covariates)
    add_primary_diagnosis!(genomicc_covariates)
    ## Process infections
    for (_, (fn, infection_name)) in severe_infections_map
        genomicc_covariates[!, infection_name] = map(fn, eachrow(genomicc_covariates))
    end
    ## Process mortality
    add_alive_at_assessment_col!(genomicc_covariates)
    ## Select cleaned columns
    DataFrames.select!(genomicc_covariates,
        :FID => :FID,
        :IID => :IID,
        :COHORT => :COHORT,
        :PRIMARY_DIAGNOSIS,
        :IS_SEVERELY_ILL,
        :AGE_YEARS_AT_RECRUITMENT => :AGE,
        :SEX_SELF_REPORTED => process_genomicc_sexes => :SEX,
        :SUPERPOPULATION,
        :AFR,
        :AMR,
        :EAS,
        :EUR,
        :SAS,
        :ALIVE_AT_ASSESSMENT,
        Symbol.(infection_names)...
    )

    # concat UKB covariates if provided
    all_covariates = ukb_covariates_file === nothing || ukb_inferred_covariates_file === nothing ?
        genomicc_covariates :
        concat_ukb_covariates(ukb_covariates_file, ukb_inferred_covariates_file, genomicc_covariates, infection_names)

    # Write to output file
    CSV.write(output_file, all_covariates, delim="\t")
    return 0
end


function make_ukb_genomicc_merge_report(;kwargs...)
    function prepend_args(script_string)
        args_string = join((string(key, " = \"", arg, "\" #hide") for (key, arg) in kwargs), "\n")
        return string(args_string, "\n", script_string)
    end
    Literate.markdown(
        joinpath(pkgdir(GenomiccWorkflows), "src", "report_ukb_genomicc_merge_template.jl"), 
        ".", 
        name="report", 
        flavor=Literate.CommonMarkFlavor(), 
        execute=true,
        preprocess=prepend_args
    )
end

function fill_chr_pvar_with_variant_id(pvar_file, variants_info_file)
    # Read variants info and pvar files
    pvar = CSV.read(pvar_file, DataFrame; delim='\t')
    variants_info = CSV.read(variants_info_file, DataFrame; delim='\t', header=["CHROM", "POS", "ID"])
    # Make mapping from position to variant ID
    variants_info_dict = Dict{Int, String}(v.POS => v.ID for v in Tables.namedtupleiterator(variants_info))
    # Update pvar with variant IDs
    pvar.ID = map(pos -> variants_info_dict[pos], pvar.POS)
    # Write updated pvar
    CSV.write(pvar_file, pvar; delim='\t', writeheader=true)

    return 0
end

function make_ukb_individuals_list(covariates_file, critical_table_file; output="ukb_eids_to_keep.txt", max_samples=nothing)
    # Read all IDs and critical IDs
    all_ids = unique(CSV.read(covariates_file, DataFrame; select=[:eid]).eid)
    critical_ids = unique(CSV.read(critical_table_file, DataFrame; select=[:eid]).eid)
    # Filter out critical IDs from all IDs and limit to max_samples if specified
    remaining_ids = setdiff(all_ids, critical_ids)
    if !isnothing(max_samples)
        remaining_ids = remaining_ids[1:max_samples]
    end
    # Write remaining eids
    open(output, "w") do f
        for id in remaining_ids
            println(f, string(id, "\t", id))
        end
    end
end