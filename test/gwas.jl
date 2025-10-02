module TestGWAS

using Test
using GenomiccWorkflows
using DataFrames
using CSV
using DelimitedFiles

PKGDIR = pkgdir(GenomiccWorkflows)
TESTDIR = joinpath(PKGDIR, "test")

function dir_contains_subdir(dir_name, subdir_name)
    subdirs = readdir(dir_name)
    if subdir_name in subdirs
        return true
    else
        first_dir = joinpath(dir_name, first(subdirs))
        if isdir(first_dir)
            return dir_contains_subdir(first_dir, subdir_name)
        else 
            return false
        end
    end
end

@testset "Test apply_filter" begin
    data = DataFrame(
        COHORT = ["GENOMICC", "GENOMICC", "UKB", "UKB", "UKB", "UKB", "UKB", "UKB", "UKB", "UKB"],
        PRIMARY_DIAGNOSIS = ["COVID-19", missing, "COVID-19", "COVID-19", "PNEUMONIA", "PNEUMONIA", "PNEUMONIA", "PNEUMONIA", "PNEUMONIA", "PNEUMONIA"],
        AGE = [25, 35, 45, 55, 65, 75, 85, missing, 50, 60]
    )
    @test GenomiccWorkflows.apply_filters(data, nothing) === data
    filter_ukb = GenomiccWorkflows.apply_filters(data, "COHORT=UKB")
    @test nrow(filter_ukb) == 8
    @test all(filter_ukb.COHORT .== "UKB")
    filter_covid_ukb = GenomiccWorkflows.apply_filters(data, "PRIMARY_DIAGNOSIS=COVID-19,COHORT=UKB")
    @test nrow(filter_covid_ukb) == 2
    @test all(filter_covid_ukb.COHORT .== "UKB")
    @test all(filter_covid_ukb.PRIMARY_DIAGNOSIS .== "COVID-19")
    @test filter_covid_ukb.AGE == [45, 55]
    filter_age_ukb = GenomiccWorkflows.apply_filters(data, "AGE>=50,AGE<=75,COHORT=UKB")
    @test nrow(filter_age_ukb) == 5
    @test all(filter_age_ukb.COHORT .== "UKB")
    @test filter_age_ukb.AGE == [55, 65, 75, 50, 60]
end

@testset "Test make-gwas-groups" begin
    tmpdir = mktempdir()
    output_prefix = joinpath(tmpdir, "gwas")
    covariates_file = joinpath(TESTDIR, "assets", "gwas", "covariates", "ukb_genomicc.covariates.csv")
    min_cases_controls = 200
    copy!(ARGS, [
        "make-gwas-groups", 
        covariates_file,
        "--groupby=SUPERPOPULATION,SEX",
        "--phenotypes=SEVERE_COVID_19",
        "--covariates=AGE,AGE_x_AGE,AGE_x_SEX,COHORT",
        "--output-prefix", output_prefix, 
        "--min-cases-controls", string(min_cases_controls)
    ])
    julia_main()

    updated_covariates = CSV.read(joinpath(tmpdir, "gwas.covariates.csv"), DataFrame)
    # Check covariate file
    expected_covariate_cols = [
        "FID", 
        "IID", 
        "AGE", 
        "SEX",
        "SUPERPOPULATION",
        "AFR",
        "SAS",
        "EAS",
        "AMR",
        "EUR",
        "COHORT",
        "SEVERE_COVID_19",
        "SEVERE_PNEUMONIA",
        "AGE_x_AGE",
        "AGE_x_SEX",
        "COHORT__GENOMICC",
        "COHORT__UKB"
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
    @test readlines(joinpath(tmpdir, "gwas.covariates_list.txt"),) == ["AGE", "AGE_x_AGE", "AGE_x_SEX", "COHORT__GENOMICC", "COHORT__UKB"]
        
    # Check groups files
    case_control_counts = sort(combine(
        groupby(updated_covariates, [:SUPERPOPULATION, :SEX, :SEVERE_COVID_19], skipmissing=true), 
        nrow), 
        :nrow
    )
    groups_failing_min_case_control_df = filter(x -> x.nrow < min_cases_controls, case_control_counts)[!, [:SUPERPOPULATION, :SEX]]
    groups_failing_min_case_control = Set(collect(zip(groups_failing_min_case_control_df.SUPERPOPULATION, groups_failing_min_case_control_df.SEX)))
    @test groups_failing_min_case_control == Set([
        ("ADMIXED", 1),
        ("ADMIXED", 0),
        ("AMR", 0),
        ("EUR", 0),
        ("AFR", 0),
        ("SAS", 0)
    ])
    for ancestry in ["AFR", "AMR", "EAS", "EUR", "SAS"]
        for sex in [0, 1]
            if (ancestry, sex) ∉ groups_failing_min_case_control
                group_key = string(ancestry, "_", sex)
                individuals = CSV.read(joinpath(tmpdir, "gwas.individuals.$group_key.SEVERE_COVID_19.txt"), DataFrame; header=["FID", "IID"])
                joined = innerjoin(updated_covariates, individuals, on = [:FID, :IID])
                @test all(==(ancestry), joined.SUPERPOPULATION)
                @test all(==(sex), joined.SEX)
                @test nrow(dropmissing(joined[!, ["SEVERE_COVID_19", "AGE", "AGE_x_AGE", "AGE_x_SEX", "COHORT__GENOMICC", "COHORT__UKB"]])) == nrow(joined)
            end
        end
    end
end

@testset "Test make-gwas-groups: no groups with filter" begin
    tmpdir = mktempdir()
    output_prefix = joinpath(tmpdir, "gwas_all")
    covariates_file = joinpath(TESTDIR, "assets", "gwas", "covariates", "ukb_genomicc.covariates.csv")
    min_cases_controls = 2500
    copy!(ARGS, [
        "make-gwas-groups", 
        covariates_file,
        "--output-prefix", output_prefix,
        "--phenotypes=SEVERE_COVID_19,SEVERE_PNEUMONIA",
        "--filters=AGE>=50,AGE<=75",
        "--covariates=AGE",
        "--min-cases-controls", string(min_cases_controls)
    ])
    julia_main()

    # Check covariate file and group 
    covariates = CSV.read(joinpath(tmpdir, "gwas_all.covariates.csv"), DataFrame)
    # SEVERE_COVID_19 is dropped because it has fewer than 2500 cases/controls
    @test !isfile(joinpath(tmpdir, "gwas_all.individuals.all.SEVERE_COVID_19.txt"))
    # The group consists in all individuals
    individuals = sort(CSV.read(
        joinpath(tmpdir, "gwas_all.individuals.all.SEVERE_PNEUMONIA.txt"), 
        DataFrame; 
        header=["FID", "IID"])
    )
    expected_individuals = sort(
        dropmissing(filter(x -> x.AGE >= 50 && x.AGE <= 75, covariates), ["SEVERE_PNEUMONIA", "AGE"]
        )[!, ["FID", "IID"]])
    @test individuals == expected_individuals

    # Check covariates list
    @test readlines(joinpath(tmpdir, "gwas_all.covariates_list.txt"),) == ["AGE"]
end

@testset "Test merge-covariates-pcs" begin
    tmpdir = mktempdir()

    covariates = DataFrame(
        FID = string.(1:10),
        IID = string.(1:10),
        SUPERPOPULATION = ["EUR", "EUR", "AFR", "AFR", "AMR", "AMR", "EAS", "EAS", "SAS", "SAS"],
        COVID_19 = [1, missing, 0, 1, 0, 1, 1, missing, 0, 1],
        AGE = rand(20:80, 10),
    )
    CSV.write(joinpath(tmpdir, "covariates.csv"), covariates, delim="\t")

    for ancestry in ["AFR", "AMR", "EAS", "EUR", "SAS"]
        for chr in 1:3
            matching_covariates = covariates[covariates.SUPERPOPULATION .== ancestry, :]
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
        "SUPERPOPULATION",
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

@testset "Test merge_chr_results" begin
    tmpdir = mktempdir()
    group = "AFR.SEVERE_COVID_19"
    chrs = 1:2
    gwas_merge_list = []
    for chr in chrs
        # Create dummy Regenie Step 2 results files
        gwas_results = DataFrame(
            CHROM = [chr],
            GENPOS = [1000],
            ID = ["chr$(chr):4132:G:A"],
            ALLELE0 = ["A"],
            ALLELE1 = ["T"],
            A1FREQ = [.5],
            N = [100],
            TEST = ["ADD"],
            BETA = [0.01],
            SE = [0.001],
            CHISQ = [10],
            LOG10P = [1.0],
            EXTRA = [""]
        )
        gwas_output_file = joinpath(tmpdir, "$group.chr$(chr).step2_SEVERE_COVID_19.regenie")
        CSV.write(gwas_output_file, gwas_results)
        push!(gwas_merge_list, gwas_output_file)
    end
    gwas_merge_list_file = joinpath(tmpdir, "gwas_merge_list.txt")
    open(gwas_merge_list_file, "w") do io
        for file in gwas_merge_list
            println(io, file)
        end
    end
    output_prefix = joinpath(tmpdir, "results.all_chr")
    copy!(ARGS, [
        "merge-chr-results",
        gwas_merge_list_file,
        "--output-prefix", output_prefix
    ])
    julia_main()
    gwas_results = CSV.read(output_prefix * ".tsv", DataFrame; delim="\t")
    expected_cols = [
        "CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", 
        "A1FREQ", "N", "TEST", "BETA", "SE", 
        "CHISQ", "LOG10P", "EXTRA"
    ]
    @test nrow(gwas_results) == 2
    @test Set(gwas_results.CHROM) == Set([1, 2])
    @test names(gwas_results) == expected_cols
end

end

true