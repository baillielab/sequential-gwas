module TestGWAS

using Test
using GenomiccWorkflows
using DataFrames
using CSV
using DelimitedFiles

PKGDIR = pkgdir(GenomiccWorkflows)
TESTDIR = joinpath(PKGDIR, "test")

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
    expected_groups = Set(["AFR", "AMR", "EAS", "EUR", "SAS"])
    # Test groups and covariates
    groups_dir = joinpath(results_dir, "call-make_covariates_and_groups", "execution")
    for ancestry in ["AFR", "AMR", "EAS", "EUR", "SAS"]
        individuals = CSV.read(joinpath(groups_dir, "grouped.individuals.$ancestry.txt"), DataFrame, header=["FID", "IID"])
        @test nrow(individuals) >= 100
    end
    covariates = CSV.read(joinpath(groups_dir, "grouped.covariates.csv"), DataFrame)
    @test "AGE_x_AGE" in names(covariates)
    covariate_list = readlines(joinpath(groups_dir, "grouped.covariates_list.txt"))
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
    ## One PCA per (group, chromosome) pair = 5 * 3 = 15
    ## These are ordered by group and chromosome
    pca_dir = joinpath(results_dir, "call-loco_pca")
    shard = 0
    for group in ["AFR", "AMR", "EAS", "EUR", "SAS"]
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
    regenie_step1_dir = joinpath(results_dir, "call-regenie_step1")
    for shard in 0:4
        execution_dir = joinpath(regenie_step1_dir, "shard-$shard", "execution")
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
    for shard in 0:14
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
        ("AMR", "1"), ("AMR", "2"), ("AMR", "3"),
        ("EUR", "1"), ("EUR", "2"), ("EUR", "3"),
        ("EAS", "1"), ("EAS", "2"), ("EAS", "3"),
        ("SAS", "1"), ("SAS", "2"), ("SAS", "3"),
        ("AFR", "1"), ("AFR", "2"), ("AFR", "3")
    ])

end

true
end