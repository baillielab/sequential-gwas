module TestQCFilesFromKGP

using Test
using GenomiccWorkflows
using CSV
using DataFrames

TESTDIR = joinpath(pkgdir(GenomiccWorkflows), "test")

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
        ("chr1", 2) => ("A", "C"),
        ("chr1", 3) => ("T", "A"),
    )
    # After merging with KGP, if the variant is not in KGP, 
    # it will have a missing VARIANT_ID
    row = (CHR_CODE = "chr1", BP_COORD=0)
    @test GenomiccWorkflows.get_action(row, kgp_info) == "DROP (VARIANT-NOT-IN-KGP)"
    # We only keep biallelic variants for now
    row = (CHR_CODE = "chr1", BP_COORD=1)
    @test GenomiccWorkflows.get_action(row, kgp_info) == "DROP (VARIANT-NON-BIALLELIC)"
    # If the reference allele is ALLELE_2 (major according to plink) and 
    # the alternative allele is ALLELE_2 (minor according to plink) then all is fine
    row = (CHR_CODE = "chr1", BP_COORD=2, ALLELE_1 = "C", ALLELE_2 = "A")
    @test GenomiccWorkflows.get_action(row, kgp_info) == "KEEP (MINOR-MAJOR-MATCHING-ALT-REF)"
    row = (CHR_CODE = "chr1", BP_COORD=3, ALLELE_1 = "A", ALLELE_2 = "T")
    @test GenomiccWorkflows.get_action(row, kgp_info) == "KEEP (MINOR-MAJOR-MATCHING-ALT-REF)"

    # There are two cases:
    ## (i) The variant is palindromic, then either:
    ##  - It wasn't properly flipped by GenomeStudio or earlier steps and needs to be flipped
    ##  - It has opposite allele frequency in our dataset (this should only happen when MAF ≈ 0.5 and we drop these cases)
    ## (ii) The variant is not palyndromic, it has an opposite allele frequency in our dataset (likely when MAF ≈ 0.5), we annotate it
    ## Non palindromic variant
    kgp_info = Dict(
        ("chr1", 1) => ("C", "A"),
        ("chr1", 2) => ("T", "A"),
        ("chr1", 3) => ("C", "G")
    )
    row = (CHR_CODE = "chr1", BP_COORD=1, ALLELE_1 = "C", ALLELE_2 = "A")
    @test GenomiccWorkflows.get_action(row, kgp_info) == "KEEP (MINOR-MAJOR-REVERSED-ALT-REF)"
    ## Palindromic variants
    ### Balanced frequencies
    row = (CHR_CODE = "chr1", BP_COORD=2, ALLELE_1 = "T", ALLELE_2 = "A", ALLELE_1_FREQ = 0.49)
    @test GenomiccWorkflows.get_action(row, kgp_info) == "DROP (PALINDROMIC-MINOR-MAJOR-REVERSED-ALT-REF-BALANCED)"
    ### Unbalanced frequencies, likely to be flipped
    row = (CHR_CODE = "chr1", BP_COORD=3, ALLELE_1 = "C", ALLELE_2 = "G", ALLELE_1_FREQ = 0.3)
    @test GenomiccWorkflows.get_action(row, kgp_info) == "FLIP (PALINDROMIC-MINOR-MAJOR-REVERSED-ALT-REF)"

    # Otherwise, the alleles are non-matching the reference at all, 
    # this can be either because:
    # - The variant hasn't been flipped properly and we check the complementary alleles against the KGP 
    # (note again, the major/minor alleles may not correspond even once flipped and we report it)
    # - There is a problem with the reported alleles in the dataset and we drop the variant (this will also be reported)
    # These variants will be dropped
    kgp_info = Dict(
        ("chr1", 1) => ("A", "C"),
        ("chr1", 2) => ("G", "A"),
        ("chr1", 3) => ("A", "C"),
        ("chr1", 4) => ("T", "G")
    )
    row = (CHR_CODE = "chr1", BP_COORD=1, ALLELE_1 = "T", ALLELE_2 = "A")
    @test GenomiccWorkflows.get_action(row, kgp_info) == "DROP (ALLELES-NOT-MATCHING-KGP)"
    row = (CHR_CODE = "chr1", BP_COORD=2, ALLELE_1 = "G", ALLELE_2 = "C")
    @test GenomiccWorkflows.get_action(row, kgp_info) == "DROP (ALLELES-NOT-MATCHING-KGP)"
    # These variants will be flipped
    row = (CHR_CODE = "chr1", BP_COORD=3, ALLELE_1 = "T", ALLELE_2 = "G")
    @test GenomiccWorkflows.get_action(row, kgp_info) == "FLIP (COMPLEMENT-MINOR-MAJOR-REVERSED-ALT-REF)"
    row = (CHR_CODE = "chr1", BP_COORD=4, ALLELE_1 = "C", ALLELE_2 = "A")
    @test GenomiccWorkflows.get_action(row, kgp_info) == "FLIP (COMPLEMENT)"
end

end

true