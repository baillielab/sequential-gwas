using SequentialGWAS
using Test
using Aqua

PKGDIR = pkgdir(SequentialGWAS)
TESTDIR = joinpath(PKGDIR, "test")

@testset "SequentialGWAS.jl" begin
    # Code quality
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(SequentialGWAS)
    end

    # Unit Tests
    @test include(joinpath(TESTDIR, "qc_from_kgp.jl"))
    @test include(joinpath(TESTDIR, "relatedness.jl"))
    @test include(joinpath(TESTDIR, "covariates.jl"))

    # Test Dataset Aggregation Workflow
    include(joinpath(TESTDIR, "combine_datasets_wgs.jl"))
    include(joinpath(TESTDIR, "combine_datasets_no_wgs.jl"))
    # Test Imputation Workflow
    include(joinpath(TESTDIR, "imputation.jl"))
    # Test UKB Merge Workflow
    include(joinpath(TESTDIR, "ukbmerge.jl"))
    # Test GWAS Workflow
    include(joinpath(TESTDIR, "gwas.jl"))
end

