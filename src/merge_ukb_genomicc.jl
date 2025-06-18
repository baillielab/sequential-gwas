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

function align_ukb_variant_ids_with_kgp_and_keep_unrelated(ukb_bed_prefix, kgp_bed_prefix; out_prefix="ukb_unrelated", threshold=0.02)
    tmpdir = mktempdir()
    # Load KGP Info
    run(`plink2 --bfile $kgp_bed_prefix --freq counts --out $kgp_bed_prefix`)
    kgp_info = SequentialGWAS.get_kgp_ref_alt(kgp_bed_prefix)
    # Load UKB variant info
    run(`plink2 --bfile $ukb_bed_prefix --freq counts --out $ukb_bed_prefix`)
    ukb_info = SequentialGWAS.load_variants_info(ukb_bed_prefix)
    # Format the chromosome and set the action column
    SequentialGWAS.set_new_columns!(ukb_info, kgp_info; threshold=threshold)
    # Make sure no variant needs to be flipped: I think (hope) UKB data is already properly aligned to the + strand
    @assert isempty(findall(startswith("FLIP"), ukb_info.ACTION)) "Some variants need to be flipped in UKB data, not sure if this is expected."
    # Write variants to keep
    variants_to_keep = ukb_info.NEW_VARIANT_ID[findall(startswith("KEEP"), ukb_info.ACTION)]
    variant_to_keep_file = joinpath(tmpdir, "variants_to_keep.txt")
    CSV.write(
        variant_to_keep_file, 
        DataFrame(VARIANT_ID=variants_to_keep),
        header=false
    )
    # Write new bim file, replace missing by chr:bp:all1:all2, they will be dropped downstream
    new_bim = DataFrames.select(ukb_info, :CHR_CODE, :NEW_VARIANT_ID, :POSITION, :BP_COORD, :ALLELE_1, :ALLELE_2)
    new_bim.CHR_CODE = replace.(new_bim.CHR_CODE, "chr" => "") # KING does not like `chr` 
    new_bim_file = joinpath(tmpdir, "new.bim")
    CSV.write(
        new_bim_file,
        new_bim, 
        header=false, 
        delim="\t"
    )
    # Drop related individual using king
    run(`king --cpus $(nthreads()) -b $ukb_bed_prefix.bed --bim $new_bim_file --fam $ukb_bed_prefix.fam --unrelated --degree 2`)
    # Drop variants using plink2
    run(`plink2 --bed $ukb_bed_prefix.bed \
    --bim $new_bim_file \
    --fam $ukb_bed_prefix.fam \
    --extract $variant_to_keep_file \
    --keep kingunrelated.txt \
    --output-chr chr26 \
    --make-bed \
    --out $out_prefix`
    )

end