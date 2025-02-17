module TestArrayGenotypesmerging

using Test
using SequentialGWAS
using DataFrames
using CSV

PKGDIR = pkgdir(SequentialGWAS)
TESTDIR = joinpath(PKGDIR, "test")
RESULTS_DIR = joinpath(PKGDIR, "results")

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
    
    # Check 1000 GP
    kgp_dir = joinpath(RESULTS_DIR, "kgp", "merged_unrelated")
    bim = SequentialGWAS.read_bim(joinpath(kgp_dir, "genotypes.merged.unrelated.bim"))
    @test unique(bim.CHR_CODE) == [string("chr", k) for k in 1:22]
    @test nrow(bim) > 100
    fam = SequentialGWAS.read_fam(joinpath(kgp_dir, "genotypes.merged.unrelated.fam"))
    unrelated = CSV.read(joinpath(kgp_dir, "kgp_unrelated_individuals.txt"), DataFrame, header=[:IID])
    @test issubset(fam.IID, unrelated.IID)

    # Check liftover
    release_2021_2023_liftover_report = CSV.read(
        joinpath(RESULTS_DIR, "array_genotypes", "lifted_over", "mock.release_2021_2023.liftover_report.csv"), 
        DataFrame
    )
    @test all(release_2021_2023_liftover_report.LIFTEDOVER .== true)
    release_r8_liftover_report = CSV.read(
        joinpath(RESULTS_DIR, "array_genotypes", "lifted_over", "mock.release_r8.liftover_report.csv"), 
        DataFrame
    )
    @test all(release_r8_liftover_report.LIFTEDOVER .== true)
    # Check qced arrays
    qced_dir = joinpath(RESULTS_DIR, "array_genotypes", "qced")
    ## Check filtered samples
    filtered_sample_files = filter(endswith("filtered_samples.csv"), readdir(qced_dir))
    map(filtered_sample_files)  do file
        filtered_samples = CSV.read(joinpath(qced_dir, file), DataFrame)
        @test names(filtered_samples) == ["IID", "FID", "FATHER_ID", "MOTHER_ID", "SEX", "PHENOTYPE"]
        @test 0 < nrow(filtered_samples) < 300 # At least some but not too many samples should be filtered
    end
    ## Check filtered variants
    filtered_variants_files = filter(endswith("filtered_variants.csv"), readdir(qced_dir))
    map(filtered_variants_files)  do file
        filtered_variants = CSV.read(joinpath(qced_dir, file), DataFrame)
        @test names(filtered_variants) == ["VARIANT_ID", "CHR_CODE", "POSITION", "BP_COORD", "ALLELE_1", "ALLELE_2"]
        @test 1 < nrow(filtered_variants) < 50 # At least some but not all variants should be filtered
    end

    # Check array genotypes shared variants
    qced_shared_variants_dir = joinpath(RESULTS_DIR, "array-genotypes", "qced_shared_variants")
    shared_variants_list = CSV.read(
        joinpath(qced_shared_variants_dir, "variants_intersection.csv"), 
        DataFrame, 
        header=[:CHR, :BP_START, :BP_END, :VARIANT_ID]
    )
    @test nrow(shared_variants_list) > 40
    bim_files = filter(endswith("shared.bim"), readdir(qced_shared_variants_dir))
    map(bim_files) do file
        bim_file = joinpath(qced_shared_variants_dir, file)
        bim = SequentialGWAS.read_bim(bim_file[1:end-4])
        
        @test Set(bim.VARIANT_ID) == Set(shared_variants_list.VARIANT_ID)
    end

    # Check wgs shared variants
    bed_shared_variants_list = CSV.read(
        joinpath(qced_shared_variants_dir, "variants_intersection.bed"), 
        DataFrame, 
        header=[:CHR, :BP_START, :BP_END]
    )
    @test nrow(bed_shared_variants_list) == nrow(shared_variants_list)
    @test bed_shared_variants_list.BP_START == shared_variants_list.BP_START .- 1
    @test bed_shared_variants_list.BP_END == shared_variants_list.BP_END .+ 1

    # Check merged genotypes
    merged_bim = SequentialGWAS.read_bim(joinpath(RESULTS_DIR, "merged", "genotypes.merged"))
    @test sort(merged_bim.VARIANT_ID) == sort(shared_variants_list.VARIANT_ID)

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