module TestArrayGenotypesmerging

using Test
using GenomiccWorkflows
using DataFrames
using CSV
using DelimitedFiles
using GenomiccUtils

PKGDIR = pkgdir(GenomiccWorkflows)
TESTDIR = joinpath(PKGDIR, "test")
RESULTS_DIR = joinpath(PKGDIR, "results_no_wgs")

@testset "Test Array Genotypes Merging" begin
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
    cmd = Cmd(["nextflow", "run", "main.nf", "-c", "test/assets/combine_datasets.no_wgs.config", "-profile", profile, "-resume"])
    cd(PKGDIR) do
        run(cmd)
    end

    # Basic checks, complete checks are implemented in `combine_datasets_wgs.jl`
    fam = read_fam(joinpath(RESULTS_DIR, "genotypes.aggregated.qced.final.fam"))
    nindividuals = nrow(fam)
    @test fam.IID == fam.FID
    @test nindividuals > 1000
    @test nrow(read_bim(joinpath(RESULTS_DIR, "genotypes.aggregated.qced.final.bim"))) > 100
    @test isfile(joinpath(RESULTS_DIR, "genotypes.aggregated.qced.final.bed"))
    covariates = CSV.read(joinpath(RESULTS_DIR, "covariates.inferred.csv"), DataFrame)
    @test Set(covariates.PLATFORM) == Set(["GSA-MD-24v3-0_A1", "GSA-MD-48v4-0_A1"])
    @test nrow(covariates) == nindividuals
    @test isfile(joinpath(RESULTS_DIR, "report.md"))
end

end
true