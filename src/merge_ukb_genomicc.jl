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

function align_ukb_variants_with_kgp_and_keep_unrelated(ukb_bed_prefix, kgp_bed_prefix; out_prefix="ukb_unrelated", threshold=0.02)
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
    run(`king --cpus $(nthreads()) -b $ukb_bed_prefix.bed --bim $new_bim_file --fam $ukb_bed_prefix.fam --unrelated --degree 2`)
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

function merge_ukb_genomicc_covariates(
    genomicc_covariates_file,
    genomicc_inferred_covariates_file,
    ukb_covariates_file,
    ukb_inferred_covariates_file,
    file_with_eids_to_exclude,;
    output_file="ukb_genomicc.covariates.csv"
    )
    # Process GenOMICC covariates
    genomicc_covariates = CSV.read(genomicc_covariates_file, DataFrame)
    genomicc_inferred_covariates = CSV.read(genomicc_inferred_covariates_file, DataFrame)
    genomicc_all_covariates = innerjoin(
        genomicc_covariates, 
        genomicc_inferred_covariates, 
        on=:genotype_file_id => :IID,
    )
    DataFrames.select!(genomicc_all_covariates,
        :genotype_file_id => :FID,
        :genotype_file_id => :IID,
        :age_years => process_genomicc_age => :AGE,
        :sex => process_genomicc_sexes => :SEX,
        :ANCESTRY_ESTIMATE,
        :AFR,
        :SAS,
        :EAS,
        :AMR,
        :EUR
    )
    genomicc_all_covariates.COHORT = fill(:GENOMICC, nrow(genomicc_all_covariates))
    # Process UKB covariates
    table_with_eids_to_exclude = CSV.read(file_with_eids_to_exclude, DataFrame, select=[:eid])
    ukb_covariates = CSV.read(ukb_covariates_file, DataFrame)
    ukb_inferred_covariates = CSV.read(ukb_inferred_covariates_file, DataFrame)
    ukb_all_covariates = innerjoin(
        ukb_covariates,
        ukb_inferred_covariates,
        on=:eid => :IID
    )
    filter!(:eid => ∉(table_with_eids_to_exclude.eid), ukb_all_covariates)
    DataFrames.select!(ukb_all_covariates,
        :eid => :FID,
        :eid => :IID,
        Symbol("34-0.0") => process_ukb_age => :AGE,
        Symbol("22001-0.0") => :SEX,
        :Superpopulation => :ANCESTRY_ESTIMATE,
        :AFR,
        :SAS,
        :EAS,
        :AMR,
        :EUR
    )
    ukb_all_covariates.COHORT = fill(:UKB, nrow(ukb_all_covariates))
    # Concatenate both datasets
    all_covariates = vcat(genomicc_all_covariates, ukb_all_covariates)
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

function make_ukb_bgen_qc_and_r2_filter_files(prefix; threshold=0.9, output=string(prefix, ".extract_list.txt"))
    variants_info = CSV.read(string(prefix, ".tsv"), DataFrame; delim='\t', header=["CHROM", "POS", "ID", "REF", "ALT", "R2"])
    pvar = CSV.read(string(prefix, ".pvar"), DataFrame; delim='\t')
    # Check REF, ALT, POS alleles match
    @assert all(variants_info.REF .== pvar.REF)
    @assert all(variants_info.ALT .== pvar.ALT)
    @assert all(variants_info.POS .== pvar.POS)
    # Write the variants_info with sufficient R2 to a file
    open(output, "w") do f
        println(filter(:R2 => >=(threshold), variants_info).ID)
        for id in filter(:R2 => >=(threshold), variants_info).ID
            println(f, id)
        end
    end
    # Update the unknown ID column in the pvar file
    pvar.ID = variants_info.ID
    tmpdir = mktempdir()
    tmpfile = joinpath(tmpdir, "temp.pvar")
    CSV.write(tmpfile, pvar; delim='\t', writeheader=true)
    mv(tmpfile, string(prefix, ".pvar"), force=true)

    return 0
end