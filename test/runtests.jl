using GenomiccWorkflows
using Test

PKGDIR = pkgdir(GenomiccWorkflows)
TESTDIR = joinpath(PKGDIR, "test")

@testset "GenOMICC Workflows" begin
    # Unit Tests
    @test include(joinpath(TESTDIR, "qc_from_kgp.jl"))
    @test include(joinpath(TESTDIR, "relatedness.jl"))
    # End-to-end Workflows
    @test include(joinpath(TESTDIR, "combine_datasets_wgs.jl"))
    @test include(joinpath(TESTDIR, "combine_datasets_no_wgs.jl"))
end

