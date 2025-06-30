# """
# 1. Extracts the variants shared between UK Biobank and the 1000 Genomes Project
# 2. Filter both datasets to include only those variants using plink
# 3. Merge the datasets using plink
# 3. Perform LD pruning using plink
# 4. Estimate ancestry using admixture
# """
# function estimate_ukb_ancestry(ukb_bed, kgp_bed)
#     ukb_bim = read_bim(joinpath(ukb_bed, ".bim"))
#     kgp_bim = read_bim(joinpath(kgp_bed, ".bim"))

#     # Step 3: Filter both datasets to shared variants
#     run(`plink --bfile $ukb_bed --extract shared.snps --make-bed --out ukb_shared`)
#     run(`plink --bfile $kgp_bed --extract shared.snps --make-bed --out kgp_shared`)
#     # Step 4: Merge datasets
#     run(`plink --bfile ukb_shared --bmerge kgp_shared.bed kgp_shared.bim kgp_shared.fam --make-bed --out $out_prefix`)

#     # Step 5: LD pruning
#     run(`plink --bfile $out_prefix --indep-pairwise 200 50 0.2 --out $out_prefix`)
#     run(`plink --bfile $out_prefix --extract $out_prefix.prune.in --make-bed --out ${out_prefix}_pruned`)

#     # Step 6: Estimate ancestry using ADMIXTURE
#     run(`admixture ${out_prefix}_pruned.bed $k`)
# end

format_chromosome!(bim) =
    bim.CHR_CODE = replace.(string.(bim.CHR_CODE), "chr" => "")

function update_variant_ids_with_map!(bim, variant_ids_map)
    unmapped_ids = Set()
    bim.VARIANT_ID = map(zip(bim.CHR_CODE, bim.BP_COORD, bim.VARIANT_ID)) do (chr, loc, old_id)
        if haskey(variant_ids_map, (chr, loc))
            variant_ids_map[(chr, loc)]
        else
            push!(unmapped_ids, old_id)
            old_id
        end
    end
    return unmapped_ids
end

function align_ukb_variant_ids_with_kgp_and_keep_unrelated(ukb_bed_prefix, kgp_bed_prefix; out_prefix="ukb_unrelated", threshold=0.02)
    tmpdir = mktempdir()
    # Load KGP variants info
    kgp_bim = SequentialGWAS.read_bim(string(kgp_bed_prefix, ".bim"))
    format_chromosome!(kgp_bim)

    kgp_variant_ids_map = Dict((chr, loc) => v_id for (chr, loc, v_id) in zip(kgp_bim.CHR_CODE, kgp_bim.BP_COORD, kgp_bim.VARIANT_ID))
    # Load UKB variants info
    ukb_bim = SequentialGWAS.read_bim(string(ukb_bed_prefix, ".bim"))
    format_chromosome!(ukb_bim)

    # Map variant IDs to KGP if possible, otherwise they will be dropped
    unmapped_ids = update_variant_ids_with_map!(ukb_bim, kgp_variant_ids_map)

    # Write variants to drop
    variant_to_drop_file = joinpath(tmpdir, "variants_to_drop.txt")
    CSV.write(
        variant_to_drop_file, 
        DataFrame(VARIANT_ID=collect(unmapped_ids)),
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
    --exclude $variant_to_drop_file \
    --keep kingunrelated.txt \
    --output-chr chr26 \
    --make-bed \
    --out $out_prefix`
    )
    return 0
end