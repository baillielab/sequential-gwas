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
    output_prefix = joinpath(tmpdir, "gwas")
    covariates_file = joinpath(TESTDIR, "assets", "gwas", "covariates", "ukb_genomicc.covariates.csv")
    min_group_size = 100
    copy!(ARGS, [
        "make-gwas-groups", 
        covariates_file,
        "--groupby=ANCESTRY_ESTIMATE,SEX",
        "--covariates=AGE,AGE_x_AGE,AGE_x_SEX",
        "--output-prefix", output_prefix, 
        "--min-group-size", string(min_group_size)
    ])
    expected_logs = [
        (:info, "Running GenOMICC Workflows: make-gwas-groups"),
        (:info, "Skipping group ADMIXED_0 because it has fewer than 100 individuals."),
        (:info, "Skipping group ADMIXED_1 because it has fewer than 100 individuals.")
    ]
    @test_logs expected_logs... julia_main()

    updated_covariates = CSV.read(joinpath(tmpdir, "gwas.covariates.csv"), DataFrame)
    # Chec newly create covariates

    # Check covariate file
    expected_covariate_cols = [
        "FID", 
        "IID", 
        "AGE", 
        "SEX",
        "ANCESTRY_ESTIMATE",
        "AFR",
        "SAS",
        "EAS",
        "AMR",
        "EUR",
        "COHORT",
        "SEVERE_COVID_19",
        "SEVERE_PNEUMONIA",
        "AGE_x_AGE",
        "AGE_x_SEX"
    ]
    @test names(updated_covariates) == expected_covariate_cols
    for row in eachrow(updated_covariates)
        @test row.AGE_x_AGE == row.AGE * row.AGE
        if row.SEX === missing
            @test row.AGE_x_SEX === missing
        else
            @test row.AGE_x_SEX == row.AGE * row.SEX
        end
    end
    # Check covariates list
    @test readlines(joinpath(tmpdir, "gwas.covariates_list.txt"),) == ["AGE", "AGE_x_AGE", "AGE_x_SEX"]
        
    # Check groups files
    for ancestry in ["EUR", "AMR", "SAS", "AFR", "EAS"]
        for sex in [0, 1]
            group_key = string(ancestry, "_", sex)
            individuals = CSV.read(joinpath(tmpdir, "gwas.individuals.$group_key.txt"), DataFrame; header=["FID", "IID"])
            joined = innerjoin(updated_covariates, individuals, on = [:FID, :IID])
            @test all(==(ancestry), joined.ANCESTRY_ESTIMATE)
            @test all(==(sex), joined.SEX)
        end
    end
end

@testset "Test make-gwas-groups: no groups" begin
    tmpdir = mktempdir()
    output_prefix = joinpath(tmpdir, "gwas_all")
    covariates_file = joinpath(TESTDIR, "assets", "gwas", "covariates", "ukb_genomicc.covariates.csv")
    min_group_size = 100
    copy!(ARGS, [
        "make-gwas-groups", 
        covariates_file,
        "--output-prefix", output_prefix, 
        "--covariates=AGE",
        "--min-group-size", string(min_group_size)
    ])
    julia_main()

    # Check covariate file and group 
    covariates = CSV.read(joinpath(tmpdir, "gwas_all.covariates.csv"), DataFrame)
    individuals = CSV.read(joinpath(tmpdir, "gwas_all.individuals.all.txt"), DataFrame; header=["FID", "IID"])
    @test covariates.FID == individuals.FID
    @test covariates.IID == individuals.IID

    # Check covariates list
    @test readlines(joinpath(tmpdir, "gwas_all.covariates_list.txt"),) == ["AGE"]
end

@testset "Test merge-covariates-pcs" begin
    tmpdir = mktempdir()

    covariates = DataFrame(
        FID = string.(1:10),
        IID = string.(1:10),
        ANCESTRY_ESTIMATE = ["EUR", "EUR", "AFR", "AFR", "AMR", "AMR", "EAS", "EAS", "SAS", "SAS"],
        COVID_19 = [1, missing, 0, 1, 0, 1, 1, missing, 0, 1],
        AGE = rand(20:80, 10),
    )
    CSV.write(joinpath(tmpdir, "covariates.csv"), covariates, delim="\t")

    for ancestry in ["AFR", "AMR", "EAS", "EUR", "SAS"]
        for chr in 1:3
            matching_covariates = covariates[covariates.ANCESTRY_ESTIMATE .== ancestry, :]
            pcs = DataFrame(
                FID = matching_covariates.FID,
                IID = matching_covariates.IID,
                PC1 = randn(2),
                PC2 = randn(2),
            )
            rename!(pcs, "FID" => "#FID")
            CSV.write(joinpath(tmpdir, "pca.$ancestry.chr$(chr)_out.eigenvec"), pcs, delim="\t")
        end
    end

    covariates_file = joinpath(tmpdir, "covariates.csv")
    pcs_prefix = joinpath(tmpdir, "pca")
    copy!(ARGS, [
        "merge-covariates-pcs", 
        covariates_file,
        pcs_prefix,
        "--output", joinpath(tmpdir, "merged_covariates_and_pcs.tsv")
    ])
    julia_main()
    merged_covariates_pcs = CSV.read(joinpath(tmpdir, "merged_covariates_and_pcs.tsv"), DataFrame)
    # Merge did not add any row
    @test nrow(merged_covariates_pcs) == 10
    # 6 new columns 2 for PCS, 3 for get_chr_out_string
    @test names(merged_covariates_pcs) == [
        "FID",
        "IID",
        "ANCESTRY_ESTIMATE",
        "COVID_19",
        "AGE",
        "CHR1_OUT_PC1",
        "CHR1_OUT_PC2",
        "CHR2_OUT_PC1",
        "CHR2_OUT_PC2",
        "CHR3_OUT_PC1",
        "CHR3_OUT_PC2"
    ]
    # missings are NA
    @test sum(merged_covariates_pcs.COVID_19 .== "NA") == 2
end

end

true