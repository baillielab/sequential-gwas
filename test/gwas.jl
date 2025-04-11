module TestGWAS

using Test
using SequentialGWAS
using DataFrames
using CSV
using DelimitedFiles

PKGDIR = pkgdir(SequentialGWAS)
TESTDIR = joinpath(PKGDIR, "test")

@testset "Test GWAS by groups" begin
    results_dir = joinpath(PKGDIR, "gwas_by_group_results")
    profile = if isinteractive()
        # The source code will be mounted in the container
        "dev"
    else
        # The source code is taken from the image
        if ENV["CI_CONTAINER"] == "docker"
            "dockerci"
        elseif ENV["CI_CONTAINER"] == "singularity"
            "singularityci"
        else
            throw(ArgumentError("Unsupported CI container"))
        end
    end
    cd(PKGDIR)
    cmd = Cmd(["nextflow", "run", "main.nf", "-entry", "GWAS", "-c", "test/assets/gwas.bygroup.config", "-profile", profile, "-resume"])
    run(cmd)
end

@testset "Test GWAS all" begin
    results_dir = joinpath(PKGDIR, "gwas_no_group_results")
    profile = if isinteractive()
        # The source code will be mounted in the container
        "dev"
    else
        # The source code is taken from the image
        if ENV["CI_CONTAINER"] == "docker"
            "dockerci"
        elseif ENV["CI_CONTAINER"] == "singularity"
            "singularityci"
        else
            throw(ArgumentError("Unsupported CI container"))
        end
    end
    cd(PKGDIR)
    cmd = Cmd(["nextflow", "run", "main.nf", "-entry", "GWAS", "-c", "test/assets/gwas.all.config", "-profile", profile, "-resume"])
    run(cmd)
end

true
end