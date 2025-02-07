module TestReadWrite

using Test
using SequentialGWAS
using CSV
using DataFrames

PKGDIR = pkgdir(SequentialGWAS)
TESTDIR = joinpath(PKGDIR, "test")

@testset "Test write_variants_intersection" begin
    tmpdir = mktempdir()
    bim_dir = joinpath(TESTDIR, "assets", "bim_files")
    copy!(
        ARGS, 
        ["write-variants-intersection", "--input-dir", bim_dir, "--output", joinpath(tmpdir, "intersection")]
    )
    julia_main()
    # Check plink file
    expected_intersection = DataFrame(
        CHR = ["chr1", "chr1", "chr1", "chr2", "chr2"],
        BP_START = [203618410, 220493812, 227748658, 10346761, 21462204],
        BP_END = [203618410, 220493812, 227748658, 10346761, 21462204],
        VARIANT_ID = ["rs12120868", "rs13376101", "GSA-rs77157578", "rs13417221", "rs7567398"]
    )
    variants_intersection = CSV.read(
        joinpath(tmpdir, "intersection.csv"), DataFrame, 
        header=[:CHR, :BP_START, :BP_END, :VARIANT_ID]
    )
    @test expected_intersection == variants_intersection
    # Check gatk file
    expected_intersection = DataFrame(
        CHR = ["chr1", "chr1", "chr1", "chr2", "chr2"],
        BP_START = [203618410, 220493812, 227748658, 10346761, 21462204],
        BP_END = [203618411, 220493813, 227748659, 10346762, 21462205],
    )
    variants_intersection = CSV.read(
        joinpath(tmpdir, "intersection.bed"), DataFrame, 
        header=[:CHR, :BP_START, :BP_END]
    )
    @test select(expected_intersection, [:CHR, :BP_START, :BP_END]) == variants_intersection
end

@test "Test complete-bim-with-ref" begin
    bim_file = joinpath(TESTDIR, "assets", "complete_bim_with_ref", "incomplete.bim")
    ref_bim_file = joinpath(TESTDIR, "assets", "complete_bim_with_ref", "variants_intersection.bim")
    copy!(
        ARGS, 
        ["complete-bim-with-ref", bim_file, ref_bim_file]
    )
    julia_main()
end
end 