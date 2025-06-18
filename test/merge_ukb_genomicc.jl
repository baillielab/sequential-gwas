module TestUKBGenOMICCMerge

using Test
using SequentialGWAS

PKGDIR = pkgdir(SequentialGWAS)
TESTDIR = joinpath(PKGDIR, "test")

# End to End Workflow run

dorun = isinteractive() || (haskey(ENV, "CI_CONTAINER") && ENV["CI_CONTAINER"] == "docker")

if dorun
    CROMWELL_PATH, CROMWELL_CONF = haskey(ENV, "CROMWELL_PATH") ? (ENV["CROMWELL_PATH"], "") : ("/Users/olabayle/cromwell/cromwell-90.jar", "-Dconfig.file=conf/cromwell.mac.conf")
    @assert isfile(CROMWELL_PATH) "Cromwell JAR file not found at $CROMWELL_PATH"

    @testset "Test Array Genotypes Merging" begin
        rc = run(`java $CROMWELL_CONF -jar $CROMWELL_PATH run wdl/ukb_merge/workflow.wdl --inputs $TESTDIR/assets/ukb_merge.json --options $TESTDIR/assets/ukb_merge_options.json`)
        @test rc.exitcode == 0
    end
end
end

true