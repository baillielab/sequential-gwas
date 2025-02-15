module TestRelatedness

using Test
using SequentialGWAS

TESTDIR = joinpath(pkgdir(SequentialGWAS), "test")

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
    kgp_unrelated_individuals(pedigree_file; outfile=outfile)
    unrelated_individuals = CSV.read(outfile, DataFrame, header=[:FID, :IID])
    @test_skip true == false
end

end

true