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
    @test include(joinpath(TESTDIR, "covariates.jl"))
    @test include(joinpath(TESTDIR, "ancestry.jl"))

    # Test Dataset Aggregation Workflow
    include(joinpath(TESTDIR, "combine_datasets_wgs.jl"))
    include(joinpath(TESTDIR, "combine_datasets_no_wgs.jl"))

    # Test Imputation Workflow
    include(joinpath(TESTDIR, "imputation.jl"))

    # Test Merging UKB and GenOMICC
    include(joinpath(TESTDIR, "merge_ukb_genomicc.jl"))

    # Test GWAS Workflow
    include(joinpath(TESTDIR, "gwas.jl"))
end

