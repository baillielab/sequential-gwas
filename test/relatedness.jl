module TestRelatedness

using Test
using GenomiccWorkflows
using CSV
using DataFrames

TESTDIR = joinpath(pkgdir(GenomiccWorkflows), "test")

@testset "Test kgp_unrelated_individuals" begin
    pedigree_file = joinpath(TESTDIR, "assets", "kgp", "20130606_g1k_3202_samples_ped_population.txt")
    tmpdir = mktempdir()
    outfile = joinpath(tmpdir, "kgp_samples_to_remove.txt")
    copy!(ARGS, [
        "get-kgp-unrelated-individuals",
        pedigree_file,
        "--output", outfile
    ])
    julia_main()
    # Check only one individual per family is kept
    pedigrees = CSV.read(pedigree_file, DataFrame)
    unrelated_individuals = CSV.read(outfile, DataFrame, header=[:IID])
    @test length(unique(pedigrees.FamilyID)) == length(unrelated_individuals.IID)
end

end

true