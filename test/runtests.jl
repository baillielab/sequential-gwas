using SequentialGWAS
using Test
using Aqua

@testset "SequentialGWAS.jl" begin
    # Code quality
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(SequentialGWAS)
    end

    # Unit Tests
    @test include("strand_flip.jl")
end

