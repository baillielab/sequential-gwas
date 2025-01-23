using Test
using SequentialGWAS

PKGDIR = pkgdir(SequentialGWAS)

@testset "Test workflow" begin
    cd(PKGDIR)
    cmd = Cmd(["nextflow", "run", "main.nf", "-c", "test/assets/workflow.config", "-profile", "ci"])
    run(cmd)
end

true