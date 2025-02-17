module TestQCFilesFromKGP

using Test
using SequentialGWAS
using CSV
using DataFrames

TESTDIR = joinpath(pkgdir(SequentialGWAS), "test")

@testset "Test get_action" begin
    # After merging with KGP, if the variant is not in KGP, 
    # it will have a missing VARIANT_ID
    row = (KGP_VARIANT_ID = missing, )
    @test SequentialGWAS.get_action(row) == "DROP (NOT-IN-KGP)"
    # We only keep biallelic variants for now
    row = (KGP_VARIANT_ID = "chr1:1:A:C", NON_BIALLELIC = true, ) 
    @test SequentialGWAS.get_action(row) == "DROP (NON-BIALLELIC)"
    # If the reference allele is ALLELE_2 (major according to plink) and 
    # the alternative allele is ALLELE_2 (minor according to plink) then all is fine
    row = (KGP_VARIANT_ID = "chr1:1:A:C", NON_BIALLELIC = false, ALLELE_1 = "C", ALLELE_2 = "A")
    @test SequentialGWAS.get_action(row) == "KEEP (MATCHING-REF-ALT)"
    row = (KGP_VARIANT_ID = "chr1:1:T:A", NON_BIALLELIC = false, ALLELE_1 = "A", ALLELE_2 = "T")
    @test SequentialGWAS.get_action(row) == "KEEP (MATCHING-REF-ALT)"

    # If the reference allele is ALLELE_1 (minor) and the alternative allele is ALLELE_1 (minor) 
    # then it means (i) if the variant is palindromic 
    # - It wasn't flipped properly by Genome Studio or previous steps
    # - It was flipped properly but has opposite major/minor alleles in our dataset as compared to the reference KGP
    # (ii) if the variant is not palindromic, it has opposite major/minor alleles in our dataset as compared to the reference KGP
    # The "opposite major/minor allele event" is unexpected, we thus flip palindromic variants and keep otherwise but will report on it
    ## Non palindromic variant
    row = (KGP_VARIANT_ID = "chr1:1:C:A", NON_BIALLELIC = false, ALLELE_1 = "C", ALLELE_2 = "A")
    @test SequentialGWAS.get_action(row) == "KEEP (REVERSE-REF-ALT)"
    ## Palindromic variants
    row = (KGP_VARIANT_ID = "chr1:1:T:A", NON_BIALLELIC = false, ALLELE_1 = "T", ALLELE_2 = "A")
    @test SequentialGWAS.get_action(row) == "FLIP (PALINDROMIC-NOT-MATCHING-KGP)"
    row = (KGP_VARIANT_ID = "chr1:1:C:G", NON_BIALLELIC = false, ALLELE_1 = "C", ALLELE_2 = "G")
    @test SequentialGWAS.get_action(row) == "FLIP (PALINDROMIC-NOT-MATCHING-KGP)"

    # Otherwise, the alleles are non-matching the reference at all, 
    # this can be either because:
    # - The variant hasn't been flipped properly and we check the complementary alleles against the KGP 
    # (note again, the major/minor alleles may not correspond even once flipped and we report it)
    # - There is a problem with the reported alleles in the dataset and we drop the variant (this will also be reported)
    # These variants will be dropped
    row = (KGP_VARIANT_ID = "chr1:1:A:C", NON_BIALLELIC = false, ALLELE_1 = "T", ALLELE_2 = "A")
    @test SequentialGWAS.get_action(row) == "DROP (ALLELES-NOT-MATCHING-KGP)"
        row = (KGP_VARIANT_ID = "chr1:1:G:A", NON_BIALLELIC = false, ALLELE_1 = "G", ALLELE_2 = "C")
    @test SequentialGWAS.get_action(row) == "DROP (ALLELES-NOT-MATCHING-KGP)"
    # These variants will be flipped
    row = (KGP_VARIANT_ID = "chr1:1:A:C", NON_BIALLELIC = false, ALLELE_1 = "T", ALLELE_2 = "G")
    @test SequentialGWAS.get_action(row) == "FLIP (COMPLEMENT-NOT-MATCHING-KGP)"
    row = (KGP_VARIANT_ID = "chr1:1:T:G", NON_BIALLELIC = false, ALLELE_1 = "C", ALLELE_2 = "A")
    @test SequentialGWAS.get_action(row) == "FLIP (COMPLEMENT)"
end

@testset "Test qc-from-kgp" begin
    outdir = mktempdir()
    release_r8_bim_file = joinpath(TESTDIR, "assets", "qc_from_kgp", "release_r8.bim")
    release_2021_2023_bim_file = joinpath(TESTDIR, "assets", "qc_from_kgp", "release_2021_2023.bim")
    release_2024_now_bim_file = joinpath(TESTDIR, "assets", "qc_from_kgp", "release_2024_now.bim")
    kgp_bim_file = joinpath(TESTDIR, "assets", "qc_from_kgp", "kgp.bim")
    copy!(ARGS, [
        "qc-from-kgp",
        "--release-r8-bim", release_r8_bim_file,
        "--release-2021-2023-bim", release_2021_2023_bim_file,
        "--release-2024-now-bim", release_2024_now_bim_file,
        "--kgp-bim", kgp_bim_file,
        "--outdir", outdir
    ])
    julia_main()
    summaries = [
        CSV.read(joinpath(outdir, "release_r8.summary.csv"), DataFrame),
        CSV.read(joinpath(outdir, "release_2021_2023.summary.csv"), DataFrame),
        CSV.read(joinpath(outdir, "release_2024_now.summary.csv"), DataFrame)
    ]
    flips = [
        CSV.read(joinpath(outdir, "release_r8.flip.txt"), DataFrame, header=[:VARIANT_ID], delim="\t"),
        CSV.read(joinpath(outdir, "release_2021_2023.flip.txt"), DataFrame, header=[:VARIANT_ID], delim="\t"),
        CSV.read(joinpath(outdir, "release_2024_now.flip.txt"), DataFrame, header=[:VARIANT_ID], delim="\t")
    ]
    new_bims = [
        SequentialGWAS.read_bim(joinpath(outdir, "release_r8.new.bim")),
        SequentialGWAS.read_bim(joinpath(outdir, "release_2021_2023.new.bim")),
        SequentialGWAS.read_bim(joinpath(outdir, "release_2024_now.new.bim"))
    ]
    old_bims = [
        SequentialGWAS.read_bim(release_r8_bim_file),
        SequentialGWAS.read_bim(release_2021_2023_bim_file),
        SequentialGWAS.read_bim(release_2024_now_bim_file)
    ]
    variants_intersection = CSV.read(joinpath(outdir, "variants_intersection.txt"), DataFrame, header=[:ID], delim="\t")
    @test length(variants_intersection.ID) > 10
    @test length(summaries) == length(flips) == length(new_bims) == 3

    for (summary, flipped, new_bim, old_bim) in zip(summaries, flips, new_bims, old_bims)
        # Dropped variants are not in the shared variants
        dropped = filter(:ACTION => startswith("DROP"), summary).VARIANT_ID
        @test length(intersect(variants_intersection.ID, dropped)) == 0
        # All flipped variants are in the shared variants
        @test all(x in variants_intersection.ID for x in flipped.VARIANT_ID)
        # Bim file has:
        # - A new VARIANT_ID column respecting the chr:pos:ref:alt format
        # - All variants from the old bim file (otherwise plink will error downstream)
        @test all(x == 4 for x in length.(split.(new_bim.VARIANT_ID, ":")))
        @test size(old_bim) == size(new_bim)
    end
end

end

true