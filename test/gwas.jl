module TestGWAS

using Test
using GenomiccWorkflows
using DataFrames
using CSV
using DelimitedFiles

PKGDIR = pkgdir(GenomiccWorkflows)
TESTDIR = joinpath(PKGDIR, "test")

@testset "Test harmonize" begin
    # If the LOG10P contains NA it will be read as a string column
    ## The NaNs are filtered for plotting
    results = DataFrame(
        LOG10P = ["NA", "1", "1", "2"],
        CHROM = [1, 1, 2, 2],
        GENPOS = [1000, 2000, 3000, 4000],
        ID = ["rs1", "rs2", "rs3", "rs4"]
    )
    harmonized_resulst = GenomiccWorkflows.harmonize(results)
    @test harmonized_resulst == DataFrame(
        CHR = ["1", "2", "2"],
        BP = [2000, 3000, 4000],
        SNP = ["rs2", "rs3", "rs4"],
        P = [0.1, 0.1, 0.01]
    )
    # If the LOG10P has no NA it will be read as a float column
    results = DataFrame(
        LOG10P = [1, 1, 2, 2],
        CHROM = [1, 1, 2, 2],
        GENPOS = [1000, 2000, 3000, 4000],
        ID = ["rs1", "rs2", "rs3", "rs4"]
    )
    harmonized_resulst = GenomiccWorkflows.harmonize(results)
    @test harmonized_resulst == DataFrame(
        CHR = ["1", "1", "2", "2"],
        BP = [1000, 2000, 3000, 4000],
        SNP = ["rs1", "rs2", "rs3", "rs4"],
        P = [0.1, 0.1, 0.01, 0.01]
    )
end


@testset "Test make-gwas-groups" begin
    tmpdir = mktempdir()
    output_prefix = joinpath(tmpdir, "gwas")
    covariates_file = joinpath(TESTDIR, "assets", "gwas", "covariates", "ukb_genomicc.covariates.csv")
    min_group_size = 100
    copy!(ARGS, [
        "make-gwas-groups", 
        covariates_file,
        "--groupby=SUPERPOPULATION,SEX",
        "--covariates=AGE,AGE_x_AGE,AGE_x_SEX,COHORT",
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
    for ancestry in ["EUR", "AMR", "SAS", "AFR", "EAS"]
        for sex in [0, 1]
            group_key = string(ancestry, "_", sex)
            individuals = CSV.read(joinpath(tmpdir, "gwas.individuals.$group_key.txt"), DataFrame; header=["FID", "IID"])
            joined = innerjoin(updated_covariates, individuals, on = [:FID, :IID])
            @test all(==(ancestry), joined.SUPERPOPULATION)
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

@testset "Test merge_regenie_chr_results" begin
    tmpdir = mktempdir()
    groups = ["AFR", "EUR"]
    chrs = 1:2
    phenotypes = ["SEVERE_COVID_19", "SEVERE_PNEUMONIA"]
    merge_list = []
    for group in groups
        for chr in chrs
            for phenotype in phenotypes
                # Create dummy Regenie Step 1 results files
                results = DataFrame(
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
                output_file = joinpath(tmpdir, "$group.chr$(chr).step2_$phenotype.regenie")
                CSV.write(output_file, results)
                push!(merge_list, output_file)
            end
        end
    end
    merge_list_file = joinpath(tmpdir, "merge_list.txt")
    open(merge_list_file, "w") do io
        for file in merge_list
            println(io, file)
        end
    end
    results_prefix = joinpath(tmpdir, "regenie.results")
    copy!(ARGS, [
        "merge-regenie-chr-results",
        merge_list_file,
        "--output-prefix", results_prefix
    ])
    julia_main()
    expected_cols = [
        "CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", 
        "A1FREQ", "N", "TEST", "BETA", "SE", 
        "CHISQ", "LOG10P", "EXTRA"
    ]
    for group in groups
        for phenotype in phenotypes
            output_file = joinpath(tmpdir, "$results_prefix.$group.$phenotype.tsv")
            results = CSV.read(output_file, DataFrame; delim="\t")
            @test nrow(results) == 2
            @test Set(results.CHROM) == Set([1, 2])
            @test names(results) == expected_cols
        end
    end
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
    expected_groups = Set(["AMR", "EAS"])
    # Test groups and covariates: only AMR and EAS have more than 2400 individuals
    groups_dir = joinpath(results_dir, "call-make_covariates_and_groups", "execution")
    for ancestry in ["AMR", "EAS"]
        individuals = CSV.read(joinpath(groups_dir, "grouped.individuals.$ancestry.txt"), DataFrame, header=["FID", "IID"])
        @test nrow(individuals) >= 2400
    end
    covariates = CSV.read(joinpath(groups_dir, "grouped.covariates.csv"), DataFrame)
    @test "AGE_x_AGE" in names(covariates)
    covariate_list = readlines(joinpath(groups_dir, "grouped.covariates_list.txt"))
    @test Set(covariate_list) == Set(["AGE", "SEX", "AGE_x_AGE"])

    # Test BED groups qced
    bed_dir = joinpath(results_dir, "call-make_group_bed_qced")
    groups = Set([])
    for shard in [0, 1]
        execution_dir = joinpath(bed_dir, "shard-$shard", "execution")
        files = readdir(execution_dir)
        fam_file = files[findfirst(endswith(".fam"), files)]
        push!(groups, splitext(fam_file)[1])
    end
    @test groups == expected_groups

    # Test LD pruning
    ld_prune_dir = joinpath(results_dir, "call-groups_ld_prune")
    groups = Set([])
    for shard in [0, 1]
        execution_dir = joinpath(ld_prune_dir, "shard-$shard", "execution")
        files = readdir(execution_dir)
        fam_file = files[findfirst(endswith(".fam"), files)]
        push!(groups, splitext(splitext(fam_file)[1])[1])
    end
    @test groups == expected_groups

    # Test LOCO PCA
    ## One PCA per (group, chromosome) pair = 5 * 3 = 15
    ## These are ordered by group and chromosome
    pca_dir = joinpath(results_dir, "call-loco_pca")
    local shard = 0
    shard = 0
    for group in ["AMR", "EAS"]
        for chr in 1:3
            execution_dir = joinpath(pca_dir, "shard-$shard", "execution")
            eigenvec_file = joinpath(execution_dir, "pca.$group.chr$(chr)_out.eigenvec")
            @test isfile(eigenvec_file)
            shard += 1
        end
    end

    # Test merge covariates and PCs
    covariates_and_pcs = CSV.read(
        joinpath(results_dir, "call-merge_covariates_and_pcs", "execution", "merged_covariates_and_pcs.tsv"), 
        DataFrame
    )
    for chr in 1:3
        for pc in 1:10
            @test "CHR$(chr)_OUT_PC$(pc)" in names(covariates_and_pcs)
        end
    end

    # Test Regenie Step 1
    regenie_step_1_dir = joinpath(results_dir, "call-regenie_step_1")
    for shard in 0:1
        execution_dir = joinpath(regenie_step_1_dir, "shard-$shard", "execution")
        files = readdir(execution_dir)
        pred_list = files[findfirst(endswith(".step1_pred.listrelative"), files)]
        phenotype_pred_files = readlines(joinpath(execution_dir, pred_list))
        phenotypes = getindex.(split.(phenotype_pred_files, " "), 1)
        @test Set(phenotypes) == Set(["SEVERE_PNEUMONIA", "SEVERE_COVID_19"])
    end

    # Test Regenie Step 2
    results_expected_cols = [
            "CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ", "N", "TEST", "BETA", "SE", "CHISQ", "LOG10P", "EXTRA"
        ]
    regenie_step2_dir = joinpath(results_dir, "call-regenie_step_2")
    ancestries_and_chrs = Set{Tuple{String, String}}([])
    for shard in 0:5
        execution_dir = joinpath(regenie_step2_dir, "shard-$shard", "execution")
        files = readdir(execution_dir)
        ## Covid-19
        covid_results_file = only(filter(f -> endswith(f, "step2_SEVERE_COVID_19.regenie"), files))
        covid_results = CSV.read(joinpath(execution_dir, covid_results_file), DataFrame)
        @test names(covid_results) == results_expected_cols
        @test nrow(covid_results) > 0
        ## Pneumonia
        pneumonia_results_file = only(filter(f -> endswith(f, "step2_SEVERE_PNEUMONIA.regenie"), files))
        pneumonia_results = CSV.read(joinpath(execution_dir, pneumonia_results_file), DataFrame)
        @test names(pneumonia_results) == results_expected_cols
        @test nrow(pneumonia_results) > 0
        # ancestry, chr 
        ancestry, chr, _ = split(covid_results_file, ".")
        push!(ancestries_and_chrs, (ancestry, chr))
    end
    @test ancestries_and_chrs == Set([
        ("AMR", "chr1"),
        ("EAS", "chr1"),
        ("AMR", "chr2"),
        ("EAS", "chr2"),
        ("AMR", "chr3"),
        ("EAS", "chr3")
    ])

    # Test merge Regenie results
    merged_results_dir = joinpath(results_dir, "call-merge_regenie_chr_results", "execution")
    for group in ["AMR", "EAS"]
        for phenotype in ["SEVERE_COVID_19", "SEVERE_PNEUMONIA"]
            output_file = joinpath(merged_results_dir, "regenie.results.$group.$phenotype.tsv")
            results = CSV.read(output_file, DataFrame; delim="\t")
            @test names(results) == results_expected_cols
            @test nrow(results) > 0
        end
    end

    # Test plots
    plots_dir = joinpath(results_dir, "call-gwas_plots")
    for shard in 0:3
        execution_dir = joinpath(plots_dir, "shard-$shard", "execution")
        files = readdir(execution_dir)
        manhattan_plot = only(filter(f -> endswith(f, ".manhattan.png"), files))
        @test isfile(joinpath(execution_dir, manhattan_plot))
        qq_file = only(filter(f -> endswith(f, ".qq.png"), files))
        @test isfile(joinpath(execution_dir, qq_file))
    end

end

end

true