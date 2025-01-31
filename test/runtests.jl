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
    include(joinpath(TESTDIR, "array_genotypes_merging.jl"))

    # Unit Tests
    include(joinpath(TESTDIR, "read_write.jl"))
end

