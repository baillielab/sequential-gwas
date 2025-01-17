using Test
using SequantialGWAS

PKGDIR = pkgdir(SequantialGWAS)

@testset "Test workflow" begin
    cd(PKGDIR)
    cmd = Cmd(["nextflow", "run", "main.nf", "-c", "test/assets/workflow.config", "-profile", "ci"])
    run(cmd)
end

true