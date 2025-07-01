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
    kgp_bim = SequentialGWAS.read_bim(string(kgp_bed_prefix, ".bim"))
    format_chromosome!(kgp_bim)

    # Create a map of variant IDs from KGP
    kgp_variant_ids_map = Dict((chr, loc) => (v_id, Set([all_1, all_2])) for (chr, loc, v_id, all_1, all_2) in 
        zip(kgp_bim.CHR_CODE, kgp_bim.BP_COORD, kgp_bim.VARIANT_ID, kgp_bim.ALLELE_1, kgp_bim.ALLELE_2))

    # Load UKB variants info
    ukb_bim = SequentialGWAS.read_bim(string(ukb_bed_prefix, ".bim"))
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

function read_covariate_and_lowercase_id(filepath; id_column=:IID)
    covariates = CSV.read(filepath, DataFrame)
    covariates[!, id_column] = lowercase.(covariates[!, id_column])
    return covariates
end

function merge_ukb_genomicc_covariates(
    genomicc_covariates,
    genomicc_inferred_covariates,
    ukb_covariates,
    ukb_inferred_covariates;
    output_file="ukb_genomicc.covariates.csv"
    )

    genomicc_covariates = "assets/rap/genomicc/covariates/a015_covariates.csv"
    genomicc_inferred_covariates = "assets/rap/genomicc/covariates/inferred_covariates.csv"

    ukb_covariates = "assets/rap/ukb/covariates/ukb_covariates.csv"
    ukb_inferred_covariates = ""

    # Process GenOMICC covariates
    genomicc_covariates = CSV.read(genomicc_covariates, DataFrame, select=[:genotype_file_id, :age_years, :sex])
    covariates[!, genotype_file_id] = lowercase.(covariates[!, genotype_file_id])

    genomicc_covariates = read_covariate_and_lowercase_id(genomicc_covariates; id_column=:genotype_file_id) 
    CSV.read(genomicc_covariates, DataFrame)

    genomicc_inferred_covariates = read_covariate_and_lowercase_id(genomicc_inferred_covariates, id_column=:IID)
    genomicc_all_covariates = innerjoin(
        genomicc_covariates, 
        genomicc_inferred_covariates, 
        on=:genotype_file_id => :IID,
    )
    

end
