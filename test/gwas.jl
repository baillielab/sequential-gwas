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
        ID = ["rs1", "rs2", "rs3", "rs4"],
        A1FREQ = ["0.1", "0.2", "NA", "0.4"]
    )
    harmonized_resulst = GenomiccWorkflows.harmonize(results)
    @test harmonized_resulst == DataFrame(
        CHR = ["1", "2"],
        BP = [2000, 4000],
        SNP = ["rs2", "rs4"],
        P = [0.1, 0.01]
    )
    # If the LOG10P has no NA it will be read as a float column
    results = DataFrame(
        LOG10P = [1, 1, 2, 2],
        CHROM = [1, 1, 2, 2],
        GENPOS = [1000, 2000, 3000, 4000],
        ID = ["rs1", "rs2", "rs3", "rs4"],
        A1FREQ = [0.1, 0.2, 0.3, 0.4]
    )
    harmonized_resulst = GenomiccWorkflows.harmonize(results, maf=0.25)
    @test harmonized_resulst == DataFrame(
        CHR = ["2", "2"],
        BP = [3000, 4000],
        SNP = ["rs3", "rs4"],
        P = [0.01, 0.01]
    )
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
                individuals = CSV.read(joinpath(tmpdir, "gwas.individuals.$group_key.txt"), DataFrame; header=["FID", "IID"])
                joined = innerjoin(updated_covariates, individuals, on = [:FID, :IID])
                @test all(==(ancestry), joined.SUPERPOPULATION)
                @test all(==(sex), joined.SEX)
                @test readlines(joinpath(tmpdir, "gwas.phenotypes.$group_key.txt")) == ["SEVERE_COVID_19"]
            end
        end
    end
end

@testset "Test make-gwas-groups: no groups" begin
    tmpdir = mktempdir()
    output_prefix = joinpath(tmpdir, "gwas_all")
    covariates_file = joinpath(TESTDIR, "assets", "gwas", "covariates", "ukb_genomicc.covariates.csv")
    min_cases_controls = 3500
    copy!(ARGS, [
        "make-gwas-groups", 
        covariates_file,
        "--output-prefix", output_prefix,
        "--phenotypes=SEVERE_COVID_19,SEVERE_PNEUMONIA",
        "--covariates=AGE",
        "--min-cases-controls", string(min_cases_controls)
    ])
    julia_main()

    # Check covariate file and group 
    covariates = CSV.read(joinpath(tmpdir, "gwas_all.covariates.csv"), DataFrame)
    # SEVERE_COVID_19 is dropped because it has fewer than 3500 cases/controls
    @test maximum(combine(groupby(covariates, :SEVERE_COVID_19, skipmissing=true), nrow).nrow) < min_cases_controls
    @test ["SEVERE_PNEUMONIA"] == readlines(joinpath(tmpdir, "gwas_all.phenotypes.all.txt"))
    # The group consists in all individuals
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
    group = "AFR"
    chrs = 1:2
    phenotype = "SEVERE_COVID_19"
    merge_list = []
    for chr in chrs
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
    output_file = joinpath(tmpdir, "$results_prefix.$group.$phenotype.tsv")
    results = CSV.read(output_file, DataFrame; delim="\t")
    expected_cols = [
        "CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", 
        "A1FREQ", "N", "TEST", "BETA", "SE", 
        "CHISQ", "LOG10P", "EXTRA"
    ]
    @test nrow(results) == 2
    @test Set(results.CHROM) == Set([1, 2])
    @test names(results) == expected_cols
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
    expected_groups = Set(["AMR", "EAS", "EUR", "AFR", "SAS"])
    # Test groups and covariates: first see which groups make the case/control constraint
    groups_prep_dir = joinpath(results_dir, "call-make_covariates_and_groups", "execution")
    covariates = CSV.read(joinpath(groups_prep_dir, "grouped.covariates.csv"), DataFrame)
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
    for group in ["EUR", "AFR", "SAS", "EAS", "AMR", "ADMIXED"]
        if group == "ADMIXED"
            @test !isfile(joinpath(groups_prep_dir, "grouped.individuals.ADMIXED.txt"))
            @test !isfile(joinpath(groups_prep_dir, "grouped.phenotypes.ADMIXED.txt"))
        else
            if group == "EAS"
                @test Set(readlines(joinpath(groups_prep_dir, "grouped.phenotypes.EAS.txt"))) == Set(["SEVERE_COVID_19", "SEVERE_PNEUMONIA"])
            else
                @test readlines(joinpath(groups_prep_dir, "grouped.phenotypes.$group.txt")) == ["SEVERE_PNEUMONIA"]
            end
            @test countlines(joinpath(groups_prep_dir, "grouped.individuals.$group.txt")) < nrow(covariates) - 100 # not all individuals
        end
    end
    covariate_list = readlines(joinpath(groups_prep_dir, "grouped.covariates_list.txt"))
    @test Set(covariate_list) == Set(["AGE", "SEX", "AGE_x_AGE"])

    # Test BED groups qced
    bed_dir = joinpath(results_dir, "call-make_group_bed_qced")
    groups = Set([])
    for shard in [0, 1, 2, 3, 4]
        execution_dir = joinpath(bed_dir, "shard-$shard", "execution")
        files = readdir(execution_dir)
        fam_file = files[findfirst(endswith(".fam"), files)]
        push!(groups, splitext(fam_file)[1])
    end
    @test groups == expected_groups

    # Test LD pruning
    ld_prune_dir = joinpath(results_dir, "call-groups_ld_prune")
    groups = Set([])
    for shard in [0, 1, 2, 3, 4]
        execution_dir = joinpath(ld_prune_dir, "shard-$shard", "execution")
        files = readdir(execution_dir)
        fam_file = files[findfirst(endswith(".fam"), files)]
        push!(groups, splitext(splitext(fam_file)[1])[1])
    end
    @test groups == expected_groups

    # Test LOCO PCA
    scatter_dirs = filter(x -> occursin("call-Scatter", x), readdir(results_dir, join=true))
    loco_pca_dir = scatter_dirs[argmin(mtime(d) for d in scatter_dirs)] # loco pca is done first
    ## One PCA per (group, chromosome) pair = 5 * 3 = 15
    ## These are ordered by group and chromosome
    pca_groups_and_chrs = Set([])
    for group_shard in [0, 1, 2, 3, 4]
        subdir = only(readdir(joinpath(loco_pca_dir, "shard-$group_shard"), join=true))
        subdir = joinpath(only(readdir(subdir, join=true)), "call-loco_pca")
        for chr_shard in [0, 1, 2]
            execution_dir = joinpath(subdir, "shard-$chr_shard", "execution")
            files = readdir(execution_dir)
            eigenvec_file = files[findfirst(endswith("eigenvec"), files)]
            _, group, chr, _ = split(eigenvec_file, ".")
            push!(pca_groups_and_chrs, (group, chr))
        end
    end
    @test pca_groups_and_chrs == Set(Iterators.product(["AFR", "AMR", "EAS", "EUR", "SAS"], ["chr1_out", "chr2_out", "chr3_out"]))

    # Test merge covariates and PCs
    covariates_and_pcs_dir = joinpath(results_dir, "call-merge_covariates_and_pcs")
    merged_covariates_groups = Set([])
    for group_shard in 0:4
        execution_dir = joinpath(covariates_and_pcs_dir, "shard-$group_shard", "execution")
        files = readdir(execution_dir)
        merged_covariates_filename = files[findfirst(endswith("merged_covariates_and_pcs.tsv"), files)]
        push!(merged_covariates_groups, split(merged_covariates_filename, ".")[1])
        covariates_and_pcs = CSV.read(joinpath(execution_dir, merged_covariates_filename), DataFrame)
        for chr in 1:3
            for pc in 1:10
                @test "CHR$(chr)_OUT_PC$(pc)" in names(covariates_and_pcs)
            end
        end
        @test "NA" in covariates_and_pcs.SEVERE_COVID_19 # missing are coded as NA
        @test "NA" in covariates_and_pcs.SEVERE_PNEUMONIA # missing are coded as NA
    end
    @test merged_covariates_groups == expected_groups

    # Test GWAS per group / phenotype
    gwas_dir = scatter_dirs[argmax(mtime(d) for d in scatter_dirs)] # loco pca is done first
    gwas_groups = Set([])
    results_expected_cols = [
            "CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ", "N", "TEST", "BETA", "SE", "CHISQ", "LOG10P", "EXTRA"
        ]
    for group_shard in 0:4
        subdir = only(readdir(joinpath(gwas_dir, "shard-$group_shard"), join=true))
        subdir = only(readdir(subdir, join=true))
        regenie_step_1_dir = joinpath(subdir, "call-regenie_step_1")
        task_dirs = readdir(subdir, join=true)
        regenie_step_2_dir = task_dirs[findfirst(x -> occursin("call-Scatter", x), task_dirs)]
        n_phenotypes = length(readdir(regenie_step_1_dir))
        observed_phenotypes = Set([])
        group = ""
        for pheno_shard in 0:n_phenotypes-1
            # Regenie step 1
            execution_dir = joinpath(regenie_step_1_dir, "shard-$pheno_shard", "execution")
            files = readdir(execution_dir)
            pred_list = files[findfirst(endswith(".step1_pred.listrelative"), files)]
            ## add group to observed gwas groups
            group = split(pred_list, ".")[1]
            push!(gwas_groups, group)
            ## add phenotype to observed phenotype for that group
            phenotype_pred_files = only(readlines(joinpath(execution_dir, pred_list)))
            phenotype = split(phenotype_pred_files, " ")[1]
            push!(observed_phenotypes, phenotype)

            # Regenie step 2 / per chromosome
            chrs_dir = joinpath(regenie_step_2_dir, "shard-$pheno_shard")
            chrs_dir = joinpath(chrs_dir, only(readdir(chrs_dir)))
            chrs_dir = joinpath(chrs_dir, only(readdir(chrs_dir)), "call-regenie_step_2")
            n_results = 0
            for chr_shard in 0:2
                execution_dir = joinpath(chrs_dir, "shard-$chr_shard", "execution")
                step_2_results = CSV.read(joinpath(execution_dir, "$group.chr$(chr_shard+1).step2_$phenotype.regenie"), DataFrame)
                @test names(step_2_results) == results_expected_cols
                @test nrow(step_2_results) > 0
                n_results += nrow(step_2_results)
            end

            # Merged results
            execution_dir = joinpath(subdir, "call-merge_regenie_chr_results", "shard-$pheno_shard", "execution")
            results = CSV.read(joinpath(execution_dir, "regenie.results.$group.$phenotype.tsv"), DataFrame; delim="\t")
            @test names(results) == results_expected_cols
            @test nrow(results) > 0
            @test n_results == nrow(results)

            # Test plots
            execution_dir = joinpath(subdir, "call-gwas_plots", "shard-$pheno_shard", "execution")
            @test isfile(joinpath(execution_dir, "gwas.plot.$group.$phenotype.manhattan.png"))
            @test isfile(joinpath(execution_dir, "gwas.plot.$group.$phenotype.qq.png"))
        end
        if group == "EAS"
            @test observed_phenotypes == Set(["SEVERE_COVID_19", "SEVERE_PNEUMONIA"])
        else
            @test observed_phenotypes == Set(["SEVERE_PNEUMONIA"])
        end
    end
    @test gwas_groups == expected_groups

end

end

true