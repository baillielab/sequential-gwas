module TestArrayGenotypesmerging

using Test
using SequentialGWAS
using DataFrames
using CSV

PKGDIR = pkgdir(SequentialGWAS)
TESTDIR = joinpath(PKGDIR, "test")
RESULTS_DIR = joinpath(PKGDIR, "results", "array-genotypes")

@testset "Test Array Genotypes Merging" begin
    profile = if isinteractive()
        # The source code will be mounted in the container
        "dev" 
    else
        # The source code is taken from the image
        "ci"
    end
    cd(PKGDIR)
    cmd = Cmd(["nextflow", "run", "main.nf", "-c", "test/assets/workflow.config", "-profile", profile, "-resume"])
    run(cmd)
end

end
true