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

    readdir(pca_dir) 
end

true
end