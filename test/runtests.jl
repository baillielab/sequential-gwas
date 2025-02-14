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

    # Array Genotypes Merging
    include(joinpath(TESTDIR, "aggregate_genetic_data.jl"))

    # Unit Tests
    @test include(joinpath(TESTDIR, "read_write.jl"))
    @test include(joinpath(TESTDIR, "qc_from_kgp.jl"))
end

