using GenomiccWorkflows
using Test
using Aqua

PKGDIR = pkgdir(GenomiccWorkflows)
TESTDIR = joinpath(PKGDIR, "test")

@testset "GenOMICC Workflows" begin
    # Code quality
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(GenomiccWorkflows)
    end

    # Unit Tests
    @test include(joinpath(TESTDIR, "qc_from_kgp.jl"))
    @test include(joinpath(TESTDIR, "relatedness.jl"))
    @test include(joinpath(TESTDIR, "ancestry.jl"))
    @test include(joinpath(TESTDIR, "fine_mapping.jl"))

    # Test Dataset Aggregation Workflow
    @test include(joinpath(TESTDIR, "combine_datasets_wgs.jl"))
    @test include(joinpath(TESTDIR, "combine_datasets_no_wgs.jl"))

    # Test Merging UKB and GenOMICC
    @test include(joinpath(TESTDIR, "merge_ukb_genomicc.jl"))

    # Test GWAS Workflow
    @test include(joinpath(TESTDIR, "gwas.jl"))
end

