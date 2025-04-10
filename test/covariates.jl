module TestCovariates

using Test
using SequentialGWAS
using DataFrames
using CSV
using YAML

PKGDIR = pkgdir(SequentialGWAS)
TESTDIR = joinpath(PKGDIR, "test")

@testset "Test make-gwas-groups" begin
    tmpdir = mktempdir()
    output_prefix = joinpath(tmpdir, "group")
    covariates_file = joinpath(TESTDIR, "assets", "genomicc", "mock.covariates.csv")
    inferred_covariates_file = joinpath(TESTDIR, "assets", "gwas", "covariates", "covariates.inferred_from_genotypes.csv")
    variables_file = joinpath(TESTDIR, "assets", "gwas", "covariates", "variables.yaml")
    min_group_size = 100
    copy!(ARGS, [
        "make-gwas-groups", 
        covariates_file,
        inferred_covariates_file,
        variables_file,
        "--output-prefix", output_prefix, 
        "--min-group-size", string(min_group_size)
    ])
    expected_logs = [
        (:info, "Running FlowOMICC: make-gwas-groups"),
        (:info, "Skipping group AFR because it has fewer than 100 individuals."),
        (:info, "Skipping group EAS because it has fewer than 100 individuals.")
    ]
    @test_logs expected_logs... julia_main()

    output_files = readdir(tmpdir)
    groups = Set(getindex.(split.(basename.(output_files), "."), 3))
    @test groups == Set(["ADMIXED", "AMR", "EUR", "SAS"])
    variables = YAML.load_file(variables_file)
    expected_covariates = ["FID", "IID", variables["covariates"]...]
    expected_phenotypes = ["FID", "IID", variables["phenotypes"]...]
    for group in groups
        covariates = CSV.read(joinpath(tmpdir, "group.covariates.$group.csv"), DataFrame)
        @test names(covariates) == expected_covariates
        @test !any(any(ismissing.(covariates[!, col])) for col in expected_covariates)
        @test nrow(covariates) >= min_group_size
        @test (covariates[!, "AGE_x_AGE"] == covariates[!, "AGE"] .* covariates[!, "AGE"])
        @test (covariates[!, "AGE_x_SEX"] == covariates[!, "AGE"] .* covariates[!, "SEX"])

        phenotypes = CSV.read(joinpath(tmpdir, "group.phenotype.$group.csv"), DataFrame)
        @test names(phenotypes) == expected_phenotypes
        @test Set(phenotypes.SEVERE_COVID_19) == Set([0, 1])

        individuals = CSV.read(joinpath(tmpdir, "group.individuals.$group.txt"), DataFrame; header=["FID", "IID"])
        @test (individuals.FID == covariates.FID == phenotypes.FID)
        @test (individuals.IID == covariates.IID == phenotypes.IID)
    end
end

@testset "Test make-gwas-groups: no groups" begin
    tmpdir = mktempdir()
    output_prefix = joinpath(tmpdir, "nogroup")
    covariates_file = joinpath(TESTDIR, "assets", "genomicc", "mock.covariates.csv")
    inferred_covariates_file = joinpath(TESTDIR, "assets", "gwas", "covariates", "covariates.inferred_from_genotypes.csv")
    variables_file = joinpath(TESTDIR, "assets", "gwas", "covariates", "variables_no_group.yaml")
    min_group_size = 100
    copy!(ARGS, [
        "make-gwas-groups", 
        covariates_file,
        inferred_covariates_file,
        variables_file,
        "--output-prefix", output_prefix, 
        "--min-group-size", string(min_group_size)
    ])
    julia_main()

    output_files = readdir(tmpdir)
    variables = YAML.load_file(variables_file)
    expected_rows = 2148
    covariates = CSV.read(joinpath(tmpdir, "nogroup.covariates.all.csv"), DataFrame)
    @test names(covariates) == ["FID", "IID", variables["covariates"]...]
    @test nrow(covariates) == expected_rows

    phenotypes = CSV.read(joinpath(tmpdir, "nogroup.phenotype.all.csv"), DataFrame)
    @test names(phenotypes) == ["FID", "IID", variables["phenotypes"]...]
    @test nrow(covariates) == expected_rows

    individuals = CSV.read(joinpath(tmpdir, "nogroup.individuals.all.txt"), DataFrame; header=["FID", "IID"])
    @test (individuals.FID == covariates.FID == phenotypes.FID)
    @test (individuals.IID == covariates.IID == phenotypes.IID)

end

end

true