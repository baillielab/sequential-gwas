module TestArrayGenotypesmerging
using Test
using SequentialGWAS

PKGDIR = pkgdir(SequentialGWAS)

@testset "Test Array Genotypes Merging" begin
    profile = Sys.ARCH == "x86_64" ? "ci" : "localarm"
    cd(PKGDIR)
    cmd = Cmd(["nextflow", "run", "main.nf", "-c", "test/assets/workflow.config", "-profile", profile, "-resume"])
    run(cmd)
end

end
true