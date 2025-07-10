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
        fam = GenomiccWorkflows.read_fam(joinpath(execution_dir, fam_file))
        push!(groups, splitext(bed_file)[1])
    end
    @test groups == Set(["AFR", "AMR", "EAS", "EUR", "SAS"])

    # Test 

    readdir(execution_dir) 
end

true
end