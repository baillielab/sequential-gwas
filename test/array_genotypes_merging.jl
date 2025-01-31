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
    # Check liftover
    @test [] == readlines(joinpath(RESULTS_DIR, "lifted_over", "mock.release_2021_2023.liftedOver.bed.unlifted"))
    @test [] == readlines(joinpath(RESULTS_DIR, "lifted_over", "mock.release_r8.liftedOver.bed.unlifted"))
    # Check qced arrays
    qced_dir = joinpath(RESULTS_DIR, "qced")
    ## Check filtered samples
    filtered_sample_files = filter(endswith("filtered_samples.csv"), readdir(qced_dir))
    filtered_samples = mapreduce(vcat, filtered_sample_files)  do file
        CSV.read(joinpath(qced_dir, file), DataFrame)
    end
    @test names(filtered_samples) == ["IID", "FID", "FATHER_ID", "MOTHER_ID", "SEX", "PHENOTYPE"]
    @test nrow(filtered_samples) > 10
    ## Check filtered variants
    filtered_variants_files = filter(endswith("filtered_variants.csv"), readdir(qced_dir))
    filtered_variants = mapreduce(vcat, filtered_variants_files)  do file
        CSV.read(joinpath(qced_dir, file), DataFrame)
    end
    @test names(filtered_variants) == ["VARIANT_ID", "CHR_CODE", "POSITION", "BP_COORD", "ALLELE_1", "ALLELE_2"]
    @test nrow(filtered_variants) > 10
    # Check shared variants
    qced_shared_variants_dir = joinpath(RESULTS_DIR, "qced_shared_variants")
    shared_variants_list = sort(readlines(joinpath(qced_shared_variants_dir, "variants_intersection.txt")))
    bim_files = filter(endswith("shared.bim"), readdir(qced_shared_variants_dir))
    map(bim_files) do file
        bim_file = joinpath(qced_shared_variants_dir, file)
        bim = SequentialGWAS.read_bim(bim_file[1:end-4])
        @test sort(bim.VARIANT_ID) == shared_variants_list
    end
    # Check merged genotypes
    merged_bim = SequentialGWAS.read_bim(joinpath(RESULTS_DIR, "merged", "genotypes.merged"))
    @test sort(merged_bim.VARIANT_ID) == shared_variants_list
    # Check merged QC
    ## Check filtered samples
    filtered_samples = CSV.read(joinpath(RESULTS_DIR, "merged_qced", "genotypes.merged.qced.filtered_samples.csv"), DataFrame)
    @test names(filtered_samples) == ["IID", "FID", "FATHER_ID", "MOTHER_ID", "SEX", "PHENOTYPE"]
    @test nrow(filtered_samples) > 10
    ## Check filtered variants
    filtered_samples = CSV.read(joinpath(RESULTS_DIR, "merged_qced", "genotypes.merged.qced.filtered_variants.csv"), DataFrame)
    @test names(filtered_samples) == ["VARIANT_ID", "CHR_CODE", "POSITION", "BP_COORD", "ALLELE_1", "ALLELE_2"]
    @test nrow(filtered_samples) == 0
end

end
true