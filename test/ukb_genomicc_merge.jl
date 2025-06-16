module TestUKBGenOMICCMerge

dorun = isinteractive() || (haskey(ENV, "CI_CONTAINER") && ENV["CI_CONTAINER"] != "docker")

if dorun
    CROMWELL_PATH = haskey(ENV, "CROMWELL_PATH") ? ENV["CROMWELL_PATH"] : "/Users/olabayle/cromwell/cromwell-90.jar"
    @assert isfile(CROMWELL_PATH) "Cromwell JAR file not found at $CROMWELL_PATH"

    using Test
    using SequentialGWAS
    PKGDIR = pkgdir(SequentialGWAS)
    TESTDIR = joinpath(PKGDIR, "test")

    @testset "Test Array Genotypes Merging" begin
        rc = run(`java -jar $CROMWELL_PATH run wdl/ukb_merge/workflow.wdl --inputs $TESTDIR/assets/ukb_merge.json`)
        @test rc.exitcode == 0
    end
end
end

true