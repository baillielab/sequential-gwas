module TestReadWrite

using Test
using SequentialGWAS
using CSV
using DataFrames

PKGDIR = pkgdir(SequentialGWAS)
TESTDIR = joinpath(PKGDIR, "test")

@testset "Test complete-bim-with-ref" begin
    tmpdir = mktempdir()
    bim_file = joinpath(TESTDIR, "assets", "complete_bim_with_ref", "incomplete.bim")
    ref_bim_file = joinpath(TESTDIR, "assets", "qc_from_kgp", "kgp.bim")
    output = joinpath(tmpdir, "complete.bim")
    copy!(
        ARGS, [
            "complete-bim-with-ref", 
            bim_file, 
            ref_bim_file,
            "--output", output
    ]
    )
    julia_main()
    @test_skip true == false # TODO
end

end

true