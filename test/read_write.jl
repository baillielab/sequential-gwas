module TestReadWrite

using Test
using SequentialGWAS

PKGDIR = pkgdir(SequentialGWAS)
TESTDIR = joinpath(PKGDIR, "test")

@testset "Test write_variants_intersection" begin
    tmpdir = mktempdir()
    bim_dir = joinpath(TESTDIR, "assets", "bim_files")
    copy!(
        ARGS, 
        ["write-variants-intersection", "--input-dir", bim_dir, "--output", joinpath(tmpdir, "intersection.txt")]
    )
    julia_main()
    variants_intersection = readlines(joinpath(tmpdir, "intersection.txt"))
    @test Set(variants_intersection) == Set(["rs12120868", "rs13376101", "GSA-rs77157578", "rs13417221", "rs7567398"])
end

end