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
    @test include(joinpath(TESTDIR, "read_write.jl"))
    @test include(joinpath(TESTDIR, "qc_from_kgp.jl"))
    @test include(joinpath(TESTDIR, "relatedness.jl"))

    # Test Genetic Data Aggregation Workflow
    include(joinpath(TESTDIR, "aggregate_genetic_data.jl"))
    # Test Imputation Workflow
    include(joinpath(TESTDIR, "imputation.jl"))
    # Test UKB Merge Workflow
    include(joinpath(TESTDIR, "ukbmerge.jl"))
    # Test GWAS Workflow
    include(joinpath(TESTDIR, "gwas.jl"))
end

