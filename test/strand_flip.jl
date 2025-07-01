module TestStrandFlip

using GenomiccWorkflows
using Test
using DataFrames

PKGDIR = pkgdir(GenomiccWorkflows)

@testset "Test make_snps_to_flip_list" begin
    manifest_file = joinpath(PKGDIR, "assets", "GSA-24v3-0_A1.csv")
    tmpdir = mktempdir()
    output = joinpath(tmpdir, "snps_to_flip.txt")
    GenomiccWorkflows.make_snps_to_flip_list(output, manifest_file)
    snps_to_flip = DataFrame(Name=readlines(output))
    snps_to_flip.ToFlip .= true
    manifest = GenomiccWorkflows.load_illumina_manifest_file(manifest_file)
    # Check that the SNPs to flip are marked with RefStrand equal to `-`
    joined = leftjoin(manifest, snps_to_flip, on =:Name)
    joined.ToFlip = coalesce.(joined.ToFlip, false)
    @test all(x == "-" for x in filter(:ToFlip => ==(true), joined).RefStrand)
    @test all(x == "+" for x in filter(:ToFlip => ==(false), joined).RefStrand)
end

end
true