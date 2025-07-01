module TestCovariates

using Test
using GenomiccWorkflows
using DataFrames
using CSV
using YAML

PKGDIR = pkgdir(GenomiccWorkflows)
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
        variables_file,
        "--inferred-covariates", inferred_covariates_file,
        "--output-prefix", output_prefix, 
        "--min-group-size", string(min_group_size)
    ])
    expected_logs = [
        (:info, "Running FlowOMICC: make-gwas-groups"),
        (:info, "Skipping group AFR because it has fewer than 100 individuals."),
        (:info, "Skipping group EAS because it has fewer than 100 individuals.")
    ]
    @test_logs expected_logs... julia_main()

    # Check covariate file
    expected_covariate_cols = ["FID", "IID", "AGE", "SEX", "AGE_x_SEX", "AGE_x_AGE", "PLATFORM__GSA-MD-24v3-0_A1", "PLATFORM__GSA-MD-48v4-0_A1"]
    covariates = CSV.read(joinpath(tmpdir, "group.covariates.csv"), DataFrame)
    @test names(covariates) == expected_covariate_cols
    @test !any(any(ismissing.(covariates[!, col])) for col in expected_covariate_cols)
    @test nrow(covariates) > 2000
    @test (covariates[!, "AGE_x_AGE"] == covariates[!, "AGE"] .* covariates[!, "AGE"])
    @test (covariates[!, "AGE_x_SEX"] == covariates[!, "AGE"] .* covariates[!, "SEX"])
    # Check phenotypes file
    expected_phenotypes = ["FID", "IID", "SEVERE_COVID_19"]
    phenotypes = CSV.read(joinpath(tmpdir, "group.phenotypes.csv"), DataFrame)
    @test names(phenotypes) == expected_phenotypes
    @test Set(phenotypes.SEVERE_COVID_19) == Set([0, 1])
    @test nrow(phenotypes) == nrow(covariates)
    # Check groups files
    inferred_covariates = CSV.read(inferred_covariates_file, DataFrame)
    for group in ["ADMIXED", "EUR", "AMR", "SAS"]
        individuals = CSV.read(joinpath(tmpdir, "group.individuals.$group.txt"), DataFrame; header=["FID", "IID"])
        joined = innerjoin(inferred_covariates, individuals, on = [:FID, :IID])
        @test all(==(group), joined.ANCESTRY_ESTIMATE)
    end
end

@testset "Test make-gwas-groups: no groups" begin
    tmpdir = mktempdir()
    output_prefix = joinpath(tmpdir, "nogroup")
    covariates_file = joinpath(TESTDIR, "assets", "genomicc", "mock.covariates.csv")
    variables_file = joinpath(TESTDIR, "assets", "gwas", "covariates", "variables_all.yaml")
    min_group_size = 100
    copy!(ARGS, [
        "make-gwas-groups", 
        covariates_file,
        variables_file,
        "--output-prefix", output_prefix, 
        "--min-group-size", string(min_group_size)
    ])
    julia_main()

    variables = YAML.load_file(variables_file)
    covariates = CSV.read(joinpath(tmpdir, "nogroup.covariates.csv"), DataFrame)
    @test names(covariates) == ["FID", "IID", variables["covariates"]...]
    @test nrow(covariates) > 2000

    phenotypes = CSV.read(joinpath(tmpdir, "nogroup.phenotypes.csv"), DataFrame)
    @test names(phenotypes) == ["FID", "IID", variables["phenotypes"]...]
    @test nrow(phenotypes) == nrow(covariates)

    individuals = CSV.read(joinpath(tmpdir, "nogroup.individuals.all.txt"), DataFrame; header=["FID", "IID"])
    @test (individuals.FID == covariates.FID == phenotypes.FID)
    @test (individuals.IID == covariates.IID == phenotypes.IID)
end

end

true