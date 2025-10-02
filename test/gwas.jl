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

# End to End Workflow run

dorun = isinteractive() || (haskey(ENV, "CI_CONTAINER") && ENV["CI_CONTAINER"] == "docker")

if dorun
    cmd_args = haskey(ENV, "CROMWELL_PATH") ?
        ["-jar", ENV["CROMWELL_PATH"]] :
        ["-Dconfig.file=conf/cromwell.mac.conf", "-jar", "/Users/olabayle/cromwell/cromwell-90.jar"]

    cmd = Cmd([
        "java", cmd_args...,
        "run", "rap_workflows/gwas/workflow.wdl",
        "--inputs", joinpath(TESTDIR, "assets", "gwas.bygroup.json"),
    ])
    rc = run(cmd)
    @test rc.exitcode == 0

    results_dirs = readdir("cromwell-executions/gwas", join=true)
    results_dir = results_dirs[argmax(mtime(d) for d in results_dirs)]
    expected_groups = Set([
        "AFR.SEVERE_PNEUMONIA", 
        "AMR.SEVERE_PNEUMONIA", 
        "EAS.SEVERE_COVID_19", 
        "EAS.SEVERE_PNEUMONIA", 
        "EUR.SEVERE_PNEUMONIA",
        "SAS.SEVERE_PNEUMONIA"
    ])
    # Test groups and covariates: first see which groups make the case/control constraint
    groups_prep_dir = joinpath(results_dir, "call-make_covariates_and_groups", "execution")
    covariates = CSV.read(joinpath(groups_prep_dir, "gwas.covariates.csv"), DataFrame)
    @test "AGE_x_AGE" in names(covariates)
    covid_19_groups_not_passing_cc_threshold = Set(filter(
        x -> x.nrow < 600,
        sort(combine(groupby(covariates, [:SUPERPOPULATION, :SEVERE_COVID_19], skipmissing=true), nrow), :nrow)
    ).SUPERPOPULATION)
    @test covid_19_groups_not_passing_cc_threshold == Set(["ADMIXED", "EUR", "AMR", "AFR", "SAS"])
    pneumonia_groups_not_passing_cc_threshold = Set(filter(
        x -> x.nrow < 600,
        sort(combine(groupby(covariates, [:SUPERPOPULATION, :SEVERE_PNEUMONIA], skipmissing=true), nrow), :nrow)
    ).SUPERPOPULATION)
    @test pneumonia_groups_not_passing_cc_threshold == Set(["ADMIXED"])
    # Now check the groups files
    for (ancestry, phenotype) in Iterators.product(
            ["AFR", "AMR", "EAS", "EUR", "SAS"],
            ["SEVERE_COVID_19", "SEVERE_PNEUMONIA"]
        )
        potential_sample_list = joinpath(groups_prep_dir, "gwas.individuals.$ancestry.$phenotype.txt")
        if "$ancestry.$phenotype" in expected_groups
            @test isfile(potential_sample_list)
        else
            @test !isfile(potential_sample_list)
        end
    end
    covariate_list = readlines(joinpath(groups_prep_dir, "gwas.covariates_list.txt"))
    @test Set(covariate_list) == Set(["AGE", "SEX", "AGE_x_AGE"])

    # Test BED groups qced
    bed_dir = joinpath(results_dir, "call-make_group_bed_qced")
    groups = Set([])
    for shard in [0, 1, 2, 3, 4, 5] # expected 6 shards for 6 sample lists
        execution_dir = joinpath(bed_dir, "shard-$shard", "execution")
        files = readdir(execution_dir)
        fam_file = files[findfirst(endswith(".fam"), files)]
        push!(groups, splitext(fam_file)[1])
    end
    @test groups == expected_groups

    # Test LD pruning
    ld_prune_dir = joinpath(results_dir, "call-groups_ld_prune")
    groups = Set([])
    for shard in [0, 1, 2, 3, 4, 5]
        execution_dir = joinpath(ld_prune_dir, "shard-$shard", "execution")
        files = readdir(execution_dir)
        fam_file = files[findfirst(endswith(".fam"), files)]
        push!(groups, splitext(splitext(fam_file)[1])[1])
    end
    @test groups == expected_groups

    # Test LOCO PCA
    scatter_dirs = filter(x -> occursin("call-Scatter", x), readdir(results_dir, join=true))
    loco_pca_dir = findfirst(
        dir_name ->  dir_contains_subdir(dir_name, "call-loco_pca"),
        scatter_dirs
    )
    loco_pca_dir = scatter_dirs[loco_pca_dir]
    ## One PCA per (group, chromosome) pair = 5 * 3 = 15
    ## These are ordered by group and chromosome
    pca_groups_and_chrs = Set([])
    for group_shard in [0, 1, 2, 3, 4, 5]
        subdir = only(readdir(joinpath(loco_pca_dir, "shard-$group_shard"), join=true))
        subdir = joinpath(only(readdir(subdir, join=true)), "call-loco_pca")
        for chr_shard in [0, 1, 2]
            execution_dir = joinpath(subdir, "shard-$chr_shard", "execution")
            files = readdir(execution_dir)
            eigenvec_file = files[findfirst(endswith("eigenvec"), files)]
            _, ancestry, phenotype, chr, _ = split(eigenvec_file, ".")
            push!(pca_groups_and_chrs, ("$ancestry.$phenotype", chr))
        end
    end
    @test pca_groups_and_chrs == Set(Iterators.product(expected_groups, ["chr1_out", "chr2_out", "chr3_out"]))

    # Test merge covariates and PCs
    covariates_and_pcs_dir = joinpath(results_dir, "call-merge_covariates_and_pcs")
    merged_covariates_groups = Set([])
    for group_shard in 0:5
        execution_dir = joinpath(covariates_and_pcs_dir, "shard-$group_shard", "execution")
        files = readdir(execution_dir)
        merged_covariates_filename = files[findfirst(endswith("merged_covariates_and_pcs.tsv"), files)]
        ancestry, phenotype, _ = split(merged_covariates_filename, ".")
        push!(merged_covariates_groups, "$ancestry.$phenotype")
        covariates_and_pcs = CSV.read(joinpath(execution_dir, merged_covariates_filename), DataFrame)
        for chr in 1:3
            for pc in 1:10
                @test "CHR$(chr)_OUT_PC$(pc)" in names(covariates_and_pcs)
            end
        end
        # Since the dataset is merged with PCs which are computed from individuals with no missing data, there is no missing data for the phenotype of interest
        @test "NA" ∉ covariates_and_pcs[!, phenotype] # missing are coded as NA
        @test eltype(covariates_and_pcs[!, phenotype]) == Int
    end
    @test merged_covariates_groups == expected_groups

    # Test REGENIE Step 1
    regenie_step_1_dir = joinpath(results_dir, "call-regenie_step_1")
    regenie_groups = Set([])
    for group_shard in 0:5
        execution_dir = joinpath(regenie_step_1_dir, "shard-$group_shard", "execution")
        files = readdir(execution_dir)
        pred_list = files[findfirst(endswith(".step1_pred.listrelative"), files)]
        ancestry, phenotype, _ = split(pred_list, ".")
        push!(regenie_groups, "$ancestry.$phenotype")
    end
    @test regenie_groups == expected_groups

    # Test REGENIE Step 2 / finemapping
    top_regenie_step_2_dir = findfirst(
        dir_name ->  dir_contains_subdir(dir_name, "call-regenie_step_2"),
        scatter_dirs
    )
    top_regenie_step_2_dir = scatter_dirs[top_regenie_step_2_dir]
    regenie_step_2_groups = Set([])
    results_expected_cols = [
            "CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ", "N", "TEST", "BETA", "SE", "CHISQ", "LOG10P", "EXTRA"
        ]
    for group_shard in 0:5
        subdir = only(readdir(joinpath(top_regenie_step_2_dir, "shard-$group_shard"), join=true))
        subdir = only(readdir(subdir, join=true))
        regenie_step_2_dir = joinpath(subdir, "call-regenie_step_2")
        for chr_shard in 0:2
            execution_dir = joinpath(regenie_step_2_dir, "shard-$chr_shard", "execution")
            files = readdir(execution_dir)
            step_2_results_file = files[findfirst(endswith(".regenie"), files)]
            ancestry, phenotype, chr_pheno, _ = split(step_2_results_file, ".")
            chr = split(chr_pheno, "_")[1]
            push!(regenie_step_2_groups, ("$ancestry.$phenotype", chr))
            step_2_results = CSV.read(joinpath(execution_dir, step_2_results_file), DataFrame)
            @test names(step_2_results) == results_expected_cols
            @test nrow(step_2_results) > 0
        end
        finemapping_dir = joinpath(subdir, "call-finemapping")
        for chr_shard in 0:2
            execution_dir = joinpath(finemapping_dir, "shard-$chr_shard", "execution")
            fp_files = filter(endswith(".tsv"), readdir(execution_dir))
            @test length(fp_files) == 2
        end
    end
    @test regenie_step_2_groups == Set(Iterators.product(expected_groups, ["chr1", "chr2", "chr3"]))

    # Test Merged gwas results
    merge_chr_results_dir = joinpath(results_dir, "call-merge_gwas_group_chr_results")
    merged_results_groups = Set([])
    for group_shard in 0:5
        execution_dir = joinpath(merge_chr_results_dir, "shard-$group_shard", "execution")
        files = readdir(execution_dir)
        gwas_merged_results_file = files[findfirst(endswith("gwas.tsv"), files)]
        ancestry, phenotype, _ = split(gwas_merged_results_file, ".")
        push!(merged_results_groups, "$ancestry.$phenotype")
    end
    @test merged_results_groups == expected_groups

    # Test Merged finemapping results
    merge_chr_results_dir = joinpath(results_dir, "call-merge_fp_group_chr_results")
    merged_results_groups = Set([])
    for group_shard in 0:5
        execution_dir = joinpath(merge_chr_results_dir, "shard-$group_shard", "execution")
        files = readdir(execution_dir)
        gwas_merged_results_file = files[findfirst(endswith("finemapping.tsv"), files)]
        ancestry, phenotype, _ = split(gwas_merged_results_file, ".")
        push!(merged_results_groups, "$ancestry.$phenotype")
    end
    @test merged_results_groups == expected_groups

    # Test Plots
    plots_dir = joinpath(results_dir, "call-gwas_group_plots")
    plots_groups = Set([])
    for group_shard in 0:5
        execution_dir = joinpath(plots_dir, "shard-$group_shard", "execution")
        files = readdir(execution_dir)
        plot_files = filter(endswith(".png"), files)
        @test length(plot_files) >= 2
        _, _, ancestry, phenotype, _ = split(first(plot_files), ".")
        push!(plots_groups, "$ancestry.$phenotype")
    end
    @test plots_groups == expected_groups

    # Test Meta-analysis
    meta_analysis_dir = joinpath(results_dir, "call-meta_analyse", "execution")
    meta_analysis_files = filter(endswith(".tsv"), readdir(meta_analysis_dir))
    meta_analysed_phenotypes = Set(String[])
    for file in meta_analysis_files
        meta_results = CSV.read(joinpath(meta_analysis_dir, file), DataFrame; delim="\t")
        @test length(unique(meta_results.ID)) == length(meta_results.ID)
        @test nrow(meta_results) > 0
        push!(meta_analysed_phenotypes, split(file, ".")[2])
        @test all(meta_results.LOG10P .>= 0.0)
    end
    @test meta_analysed_phenotypes == Set(["SEVERE_PNEUMONIA", "SEVERE_COVID_19"])

    # Test Meta-analysis finemapping
    meta_fp_analysis_dir = joinpath(results_dir, "call-merge_fp_meta_chr_results")
    for shard in [0, 1] # 2 shards for 2 phenotypes
        execution_dir = joinpath(meta_fp_analysis_dir, "shard-$shard", "execution")
        @test findfirst(endswith("finemapping.tsv"), readdir(execution_dir)) !== nothing
    end

    # Test Meta-analysis plots
    meta_plots_dir = joinpath(results_dir, "call-gwas_meta_plots")
    for shard in [0, 1] # 2 shards for 2 phenotypes
        execution_dir = joinpath(meta_plots_dir, "shard-$shard", "execution")
        files = readdir(execution_dir)
        plot_files = filter(endswith(".png"), files)
        @test length(plot_files) == 2
    end
end

end

true