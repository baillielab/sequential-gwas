module TestQCFilesFromKGP

using Test
using GenomiccWorkflows
using CSV
using DataFrames

TESTDIR = joinpath(pkgdir(GenomiccWorkflows), "test")

@testset "Test update_bim_with_minor_major_info!" begin
    bim = DataFrame(
        CHR_CODE = ["chr1", "chr1", "chr2"],
        VARIANT_ID = ["chr1:100:A:T", "chr1:200:G:C", "chr2:300:T:A"],
        POSITION = [100, 200, 300],
        BP_COORD = [100, 200, 300],
        ALLELE_1 = ["A", "G", "T"],
        ALLELE_2 = ["T", "C", "A"]
    )
    # order in acount may not be the same as in the bim file
    acount = DataFrame(
        ID = ["chr1:200:G:C", "chr2:300:T:A", "chr1:100:A:T"],
        REF = ["C", "A", "A"],
        ALT = ["G", "T", "T"],
        ALT_CTS = [60, 50, 10],
        OBS_CT = [100, 100, 100]
    )
    GenomiccWorkflows.update_bim_with_minor_major_info!(bim, acount)
    @test bim.MINOR_ALLELE == ["T", "C", "T"]
    @test bim.MAJOR_ALLELE == ["A", "G", "A"]
    @test bim.MINOR_ALLELE_FREQ == [0.1, 0.4, 0.5]
end

@testset "Test kgp_minor_major_alleles_by_position_from_df" begin
    kgp_info = DataFrame(
        CHR_CODE = ["chr1", "chr2", "chr2"],
        VARIANT_ID = ["chr1:100:A:T", "chr1:200:G:C", "chr2:200:T:A"],
        POSITION = [100, 200, 200],
        BP_COORD = [100, 200, 200],
        ALLELE_1 = ["A", "G", "T"],
        ALLELE_2 = ["T", "C", "A"],
        MINOR_ALLELE = ["T", "C", "T"],
        MAJOR_ALLELE = ["A", "G", "A"],
        MINOR_ALLELE_FREQ = [0.1, 0.4, 0.5]
    )
    kgp_info_dict = GenomiccWorkflows.kgp_minor_major_alleles_by_position_from_df(kgp_info)
    @test kgp_info_dict == Dict(
        ("chr2", 200) => nothing, # multi-allelic variant
        ("chr1", 100) => ("T", "A", 0.1, "chr1:100:A:T")
    )
end


@testset "Test write_release_samples_to_drop" begin
    tmpdir = mktempdir()
    prefixes_and_fams = (
        release_2024_now = (
            joinpath(tmpdir, "release_2024_now"),
            DataFrame(FID = ["F1", "F2"], IID = ["I5", "I6"])
        ),
        release_2021_2023 = (
            joinpath(tmpdir, "release_2021_2023"),
            DataFrame(FID = ["F1", "F2"], IID = ["I3", "I5"])
        ),
        release_r8 = (
            joinpath(tmpdir, "release_r8"),
            DataFrame(FID = ["F1", "F2"], IID = ["I1", "I3"])
        )
    )
    # No WGS samples
    wgs_samples_file = joinpath(tmpdir, "wgs_samples.txt")
    touch(wgs_samples_file)
    GenomiccWorkflows.write_release_samples_to_drop(prefixes_and_fams, wgs_samples_file)
    @test nrow(CSV.read(joinpath(tmpdir, "release_2024_now.samples_to_drop.txt"), DataFrame, header=[:FID, :IID])) == 0
    @test CSV.read(joinpath(tmpdir, "release_2021_2023.samples_to_drop.txt"), DataFrame, header=[:FID, :IID]) == DataFrame(FID = ["F2"], IID=["I5"])
    @test CSV.read(joinpath(tmpdir, "release_r8.samples_to_drop.txt"), DataFrame, header=[:FID, :IID]) == DataFrame(FID = ["F2"], IID=["I3"])
    # WGS samples
    write(wgs_samples_file, "I5\nI1")
    GenomiccWorkflows.write_release_samples_to_drop(prefixes_and_fams, wgs_samples_file)
    @test CSV.read(joinpath(tmpdir, "release_2024_now.samples_to_drop.txt"), DataFrame, header=[:FID, :IID]) == DataFrame(FID = ["F1"], IID=["I5"])
    @test CSV.read(joinpath(tmpdir, "release_2021_2023.samples_to_drop.txt"), DataFrame, header=[:FID, :IID]) == DataFrame(FID = ["F2"], IID=["I5"])
    @test CSV.read(joinpath(tmpdir, "release_r8.samples_to_drop.txt"), DataFrame, header=[:FID, :IID]) == DataFrame(FID = ["F1", "F2"], IID=["I1", "I3"])
end


@testset "Test get_action" begin
    kgp_info = Dict(
        ("chr1", 1) => nothing,
        ("chr1", 2) => ("A", "C", 0.1, "v_id"),
        ("chr1", 3) => ("T", "A", 0.1, "v_id"),
    )
    # After merging with KGP, if the variant is not in KGP, 
    # it will have a missing VARIANT_ID
    row = (CHR_CODE = "chr1", BP_COORD=0)
    @test GenomiccWorkflows.get_action(row, kgp_info) == "DROP (VARIANT-NOT-IN-KGP)"
    # We only keep biallelic variants for now
    row = (CHR_CODE = "chr1", BP_COORD=1)
    @test GenomiccWorkflows.get_action(row, kgp_info) == "DROP (VARIANT-NON-BIALLELIC)"
    # If minor/major alleles are concordant with KGP then all is fine
    row = (CHR_CODE = "chr1", BP_COORD=2, MINOR_ALLELE = "A", MAJOR_ALLELE = "C")
    @test GenomiccWorkflows.get_action(row, kgp_info) == "KEEP (MINOR-MAJOR-MATCHING-KGP)"
    row = (CHR_CODE = "chr1", BP_COORD=3, MINOR_ALLELE = "T", MAJOR_ALLELE = "A")
    @test GenomiccWorkflows.get_action(row, kgp_info) == "KEEP (MINOR-MAJOR-MATCHING-KGP)"

    # There are two cases:
    ## (i) The variant is palindromic, then either:
    ##  - It wasn't properly flipped by GenomeStudio or earlier steps and needs to be flipped
    ##  - It has opposite allele frequency in our dataset (this should only happen when MAF ≈ 0.5 and we drop these cases)
    ## (ii) The variant is not palyndromic, it has an opposite allele frequency in our dataset (likely when MAF ≈ 0.5), we annotate it
    ## Non palindromic variant
    kgp_info = Dict(
        ("chr1", 1) => ("C", "A", 0.1, "v_id"),
        ("chr1", 2) => ("T", "A", 0.1, "v_id"),
        ("chr1", 3) => ("C", "G", 0.1, "v_id")
    )
    row = (CHR_CODE = "chr1", BP_COORD=1, MINOR_ALLELE = "A", MAJOR_ALLELE = "C")
    @test GenomiccWorkflows.get_action(row, kgp_info) == "KEEP (MINOR-MAJOR-REVERSED-KGP)"
    ## Palindromic variants
    ### Balanced frequencies
    row = (CHR_CODE = "chr1", BP_COORD=2, MINOR_ALLELE = "A", MAJOR_ALLELE = "T", MINOR_ALLELE_FREQ = 0.49)
    @test GenomiccWorkflows.get_action(row, kgp_info) == "DROP (PALINDROMIC-MINOR-MAJOR-REVERSED-KGP-BALANCED)"
    ### Unbalanced frequencies, likely to be flipped
    row = (CHR_CODE = "chr1", BP_COORD=3, MINOR_ALLELE = "G", MAJOR_ALLELE = "C", MINOR_ALLELE_FREQ = 0.3)
    @test GenomiccWorkflows.get_action(row, kgp_info) == "FLIP (PALINDROMIC-MINOR-MAJOR-REVERSED-KGP-UNBALANCED)"

    # Otherwise, the alleles are non-matching the reference at all, 
    # this can be either because:
    # - The variant hasn't been flipped properly and we check the complementary alleles against the KGP 
    # (note again, the major/minor alleles may not correspond even once flipped and we report it)
    # - There is a problem with the reported alleles in the dataset and we drop the variant (this will also be reported)
    # These variants will be dropped
    kgp_info = Dict(
        ("chr1", 1) => ("A", "C", 0.1, "v_id"),
        ("chr1", 2) => ("G", "A", 0.1, "v_id"),
        ("chr1", 3) => ("A", "C", 0.1, "v_id"),
        ("chr1", 4) => ("T", "G", 0.1, "v_id")
    )
    row = (CHR_CODE = "chr1", BP_COORD=1, MINOR_ALLELE = "A", MAJOR_ALLELE = "T")
    @test GenomiccWorkflows.get_action(row, kgp_info) == "DROP (ALLELES-NOT-MATCHING-KGP)"
    row = (CHR_CODE = "chr1", BP_COORD=2, MINOR_ALLELE = "C", MAJOR_ALLELE = "G")
    @test GenomiccWorkflows.get_action(row, kgp_info) == "DROP (ALLELES-NOT-MATCHING-KGP)"
    # These variants will be flipped
    row = (CHR_CODE = "chr1", BP_COORD=3, MINOR_ALLELE = "G", MAJOR_ALLELE = "T")
    @test GenomiccWorkflows.get_action(row, kgp_info) == "FLIP (COMPLEMENT-KGP-MINOR-MAJOR-REVERSED)"
    row = (CHR_CODE = "chr1", BP_COORD=4, MINOR_ALLELE = "A", MAJOR_ALLELE = "C")
    @test GenomiccWorkflows.get_action(row, kgp_info) == "FLIP (COMPLEMENT-KGP)"
end

@testset "Test set_new_id_column!" begin
    kgp_info = Dict(
        ("chr1", 1) => ("A", "C", 0.1, "v_id_1"),
        ("chr1", 2) => ("T", "A", 0.1, "v_id_2"),
        ("chr1", 3) => nothing
    )
    df = DataFrame(
        CHR_CODE = ["chr1", "chr1", "chr1", "chr2"],
        BP_COORD = [1, 2, 3, 4],
        VARIANT_ID = ["old_id_1", "old_id_2", "old_id_3", "old_id_4"]
    )
    GenomiccWorkflows.set_new_id_column!(df, kgp_info)
    @test df.NEW_VARIANT_ID == ["v_id_1", "v_id_2", "old_id_3", "old_id_4"]
end

end

true