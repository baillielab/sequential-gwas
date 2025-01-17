using SequantialGWAS
using Test
using Aqua

@testset "SequantialGWAS.jl" begin
    # @testset "Code quality (Aqua.jl)" begin
    #     Aqua.test_all(SequantialGWAS)
    # end

    # Test workflow
    @test include("workflow.jl")
end

