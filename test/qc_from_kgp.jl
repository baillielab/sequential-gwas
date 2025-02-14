module TestQCFromKGP

using Test
using SequentialGWAS
using CSV
using DataFrames

TESTDIR = joinpath(pkgdir(SequentialGWAS), "test")

@testset "Test qc-from-kgp"
    outdir = mktempdir()
    release_r8_bim_file = joinpath(TESTDIR, "assets", "qc_from_kgp", "release_r8.bim")
    release_2021_2023_bim_file = joinpath(TESTDIR, "assets", "qc_from_kgp", "release_2021_2023.bim")
    release_2024_now_bim_file = joinpath(TESTDIR, "assets", "qc_from_kgp", "release_2024_now.bim")
    kgp_bim_file = joinpath(TESTDIR, "assets", "qc_from_kgp", "kgp.bim")
    copy!(ARGS, [
        "qc-from-kgp",
        "--release-r8-bim", release_r8_bim_file,
        "--release-2021-2023-bim", release_2021_2023_bim_file,
        "--release-2024-now-bim", release_2024_now_bim_file,
        "--kgp-bim", kgp_bim_file,
        "--outdir", outdir
    ])
    julia_main()
    summaries = [
        CSV.read(joinpath(outdir, "release_r8.summary.csv"), DataFrame),
        CSV.read(joinpath(outdir, "release_2021_2023.summary.csv"), DataFrame),
        CSV.read(joinpath(outdir, "release_2024_now.summary.csv"), DataFrame)
    ]
    flips = [
        CSV.read(joinpath(outdir, "release_r8.flip.txt"), DataFrame, header=false),
        CSV.read(joinpath(outdir, "release_2021_2023.flip.txt"), DataFrame, header=false),
        CSV.read(joinpath(outdir, "release_2024_now.flip.txt"), DataFrame, header=false)
    ]
    new_bims = [
        SequentialGWAS.read_bim(joinpath(outdir, "release_r8.new.bim")),
        SequentialGWAS.read_bim(joinpath(outdir, "release_2021_2023.new.bim")),
        SequentialGWAS.read_bim(joinpath(outdir, "release_2024_now.new.bim"))
    ]
    variants_intersection = CSV.read(joinpath(outdir, "variants_intersection.txt"), DataFrame, header=[:ID], delim="\t")
end

true