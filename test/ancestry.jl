module TestAncestry

using Test
using GenomiccWorkflows
using CSV
using DataFrames

TESTDIR = joinpath(pkgdir(GenomiccWorkflows), "test")

@testset "Test processing input/output of scope" begin
    # Test formatting of input file for SCOPE
    input_file = joinpath(TESTDIR, "assets", "scope_inputs", "kgp.frq.strat")
    output_file, io = mktemp()
    GenomiccWorkflows.format_stratified_freqs(input_file, output_file)
    close(io)
    scope_freqs = CSV.read(output_file, DataFrame; delim="\t")
    @test Set(scope_freqs.CLST) == Set([1, 2, 3, 4, 5])
    @test scope_freqs.MAF isa Vector{Float64}
    @test all(scope_freqs.MAC .== "N")
    @test all(scope_freqs.NCHROBS .== "N")
    # Test reading scope estimates
    scope_ancestry_file = "test/assets/scope_inputs/scope_resultsQhat.txt"
    n_indiv = 150
    fam = DataFrame(FID=1:150, IID=1:150)
    Q = GenomiccWorkflows.read_scope_ancestry_estimates(n_indiv, scope_ancestry_file)
    @test size(Q) == (n_indiv, 5)
    GenomiccWorkflows.assign_scope_ancestry_estimates!(fam, Q; threshold=0.6)
    # Check first few lines
    @test fam.Superpopulation[1:6] == ["ADMIXED", "ADMIXED", "EUR", "ADMIXED", "AFR", "AFR"]
end

end

true