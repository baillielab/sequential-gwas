using Test
using SequantialGWAS

PKGDIR = pkgdir(SequantialGWAS)

@testset "Test workflow" begin
    cmd = Cmd(["nextflow", "run", "$PKGDIR/main.nf", "-c", "$PKGDIR/test/assets/workflow.config", "-profile", "ci"])
    run(cmd)
end

true