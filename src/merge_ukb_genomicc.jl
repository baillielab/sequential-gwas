function get_severe_infections_map()
    severe_infections_map = Dict{Any, Any}(
        "covid-19" => (is_severe_covid_19, "SEVERE_COVID_19"),
    )
    for (raw_name, clean_name) in [
        "influenza virus" => "SEVERE_INFLUENZA",
        "pneumonia with radiographic changes at presentation to critical care" => "SEVERE_PNEUMONIA",
        "pancreatitis of any aetiology" => "SEVERE_PANCREATITIS",
        "rsv (respiratory syncytial virus) infection" => "SEVERE_RSV",
        "soft tissue infections causing systemic sepsis" => "SEVERE_SOFT_TISSUE_INFECTION",
        "ecls" => "SEVERE_ECLS",
        "reaction to vaccination" => "SEVERE_REACTION_TO_VACCINATION"]
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
    is_severe_covid_19(row)

A specific set of rules to process the severity of COVID-19 based on the cohort and diagnosis.
"""
function is_severe_covid_19(row)
    if row.PRIM_DIAGNOSIS_ODAP == "covid-19"
        cohort = row.COHORT
        if cohort == "genomicc_severe"
            return 1
        elseif cohort == "gen_int_pakistan"
            return 1
        elseif cohort == "genomicc_mild"
            return 0
        elseif cohort == "genomicc_react"
            return 0
        elseif cohort == "isaric4c"
            if row.ISARIC_MAX_SEVERITY_SCORE == "NA"
                return missing
            else
                score = parse(Int, row.ISARIC_MAX_SEVERITY_SCORE)
                return score >= 4 ? 1 : 0
            end
        else
            throw(ArgumentError("Unknown cohort: $cohort"))
        end
    else
        return missing
    end
end

"""
    is_severe_infection(row, infection_name)

An individual is considered severe for an infection if it had the infection and is part of the genomicc_severe cohort. 
If it does not have the infection it is ignored.
If an individual has the infection but is not part of genomicc_severe, we do not know where they come from and how to process them so we throw.
"""
function is_severe_infection(row, infection_name)
    if row.PRIM_DIAGNOSIS_ODAP == infection_name
        cohort = row.COHORT
        if cohort == "genomicc_severe" || cohort == "gen_int_pakistan"
            return 1
        else
            throw(ArgumentError("Unknown cohort: $cohort, while processing infection: $infection_name"))
        end
    else
        return missing
    end
end

function merge_ukb_genomicc_covariates(
    genomicc_covariates_file,
    ukb_covariates_file,
    ukb_inferred_covariates_file;
    output_file="ukb_genomicc.covariates.csv"
    )
    severe_infections_map = get_severe_infections_map()
    infection_names = last.(values(severe_infections_map))
    # Process GenOMICC covariates
    genomicc_covariates = CSV.read(genomicc_covariates_file, DataFrame)
    ## covid-19 values may be prefixed by cohort
    genomicc_covariates.PRIM_DIAGNOSIS_ODAP = replace.(genomicc_covariates.PRIM_DIAGNOSIS_ODAP, 
        "isaric4c covid-19" => "covid-19",
        "react covid-19" => "covid-19",
        "mild covid-19" => "covid-19",
    )
    ## Process infections
    for (_, (fn, infection_name)) in severe_infections_map
        genomicc_covariates[!, infection_name] = map(fn, eachrow(genomicc_covariates))
    end
    DataFrames.select!(genomicc_covariates,
        :FID => :FID,
        :IID => :IID,
        :COHORT => :COHORT,
        :AGE_YEARS_AT_RECRUITMENT => :AGE,
        :SEX_SELF_REPORTED => process_genomicc_sexes => :SEX,
        :SUPERPOPULATION,
        Symbol.(infection_names)...
    )
    # Process UKB covariates
    ukb_covariates = CSV.read(ukb_covariates_file, DataFrame)
    ukb_inferred_covariates = CSV.read(ukb_inferred_covariates_file, DataFrame)
    ukb_all_covariates = innerjoin(
        ukb_covariates,
        ukb_inferred_covariates,
        on=:eid => :IID
    )
    DataFrames.select!(ukb_all_covariates,
        :eid => :FID,
        :eid => :IID,
        Symbol("34-0.0") => process_ukb_age => :AGE,
        Symbol("22001-0.0") => :SEX,
        :Superpopulation => :SUPERPOPULATION,
    )
    for (_, (_, infection_name)) in severe_infections_map
        ukb_all_covariates[!, infection_name] = zeros(Int, nrow(ukb_all_covariates))
    end
    ukb_all_covariates.COHORT = fill(:ukbiobank, nrow(ukb_all_covariates))
    # Concatenate both datasets
    all_covariates = vcat(genomicc_covariates, ukb_all_covariates)
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