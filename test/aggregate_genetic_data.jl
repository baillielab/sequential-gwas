module TestArrayGenotypesmerging

using Test
using SequentialGWAS
using DataFrames
using CSV
using DelimitedFiles

PKGDIR = pkgdir(SequentialGWAS)
TESTDIR = joinpath(PKGDIR, "test")
RESULTS_DIR = joinpath(PKGDIR, "results")

@testset "Test Array Genotypes Merging" begin
    profile = if isinteractive()
        # The source code will be mounted in the container
        "dev" 
    else
        # The source code is taken from the image
        if ENV["CI_CONTAINER"] == "docker"
            "dockerci"
        elseif ENV["CI_CONTAINER"] == "singularity"
            "singularityci"
        else
            throw(ArgumentError("Unsupported CI container"))
        end
    end
    cd(PKGDIR)
    cmd = Cmd(["nextflow", "run", "main.nf", "-c", "test/assets/workflow.config", "-profile", profile, "-resume"])
    run(cmd)

    # Check 1000 GP
    kgp_dir = joinpath(RESULTS_DIR, "kgp", "merged_unrelated")
    bim = SequentialGWAS.read_bim(joinpath(kgp_dir, "kgp.merged.unrelated.bim"))
    @test unique(bim.CHR_CODE) == [string("chr", k) for k in 1:22]
    @test nrow(bim) > 100
    fam = SequentialGWAS.read_fam(joinpath(kgp_dir, "kgp.merged.unrelated.fam"))
    unrelated = CSV.read(joinpath(kgp_dir, "kgp_unrelated_individuals.txt"), DataFrame, header=[:IID])
    @test issubset(fam.IID, unrelated.IID)

    # Check liftover
    @test read(joinpath(RESULTS_DIR, "array_genotypes", "lifted_over", "mock.release_2021_2023.liftedOver.bed.unlifted")) == []
    @test read(joinpath(RESULTS_DIR, "array_genotypes", "lifted_over", "mock.release_r8.liftedOver.bed.unlifted")) == []

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
        @test 0 < nrow(filtered_variants) < 50 # At least some but not all variants should be filtered
    end

    # Check 1000 GP based QC
    kgp_qc_files_dir = joinpath(RESULTS_DIR, "array_genotypes", "qc_files_from_kgp")
    flipped_shared_dir = joinpath(RESULTS_DIR, "array_genotypes", "flipped_and_shared")
    shared_variants = readlines(joinpath(kgp_qc_files_dir, "variants_intersection.txt"))
    @test length(shared_variants) > 10
    ## New bim files
    flipped_bim_files = joinpath.(flipped_shared_dir, ["release-r8.flipped.shared.bim", "release-2021-2023.flipped.shared.bim", "release-2024-now.flipped.shared.bim"])
    flip_files = joinpath.(kgp_qc_files_dir, ["mock.release_r8.liftedOver.qced.flip.txt", "mock.release_2021_2023.liftedOver.qced.flip.txt", "mock.release_2024_now.qced.flip.txt"])
    unflipped_bim_files = joinpath.(kgp_qc_files_dir, ["mock.release_r8.liftedOver.qced.new.bim", "mock.release_2021_2023.liftedOver.qced.new.bim", "mock.release_2024_now.qced.new.bim"])

    for (flipped_bim_file, unflipped_bim_file, flip_file) in zip(flipped_bim_files, unflipped_bim_files, flip_files)
        flipped_bim = SequentialGWAS.read_bim(flipped_bim_file)
        # Check chr:pos:ref:alt format
        @test all(length.(split.(flipped_bim.VARIANT_ID, ":")) .== 4)
        # Check all shared variants are present in the bim file
        @test sort(shared_variants) == sort(flipped_bim.VARIANT_ID)
        # Check variants have been flipped
        unflipped_bim = SequentialGWAS.read_bim(unflipped_bim_file)
        select!(unflipped_bim, 
            :VARIANT_ID, 
            :ALLELE_1 => :UNFLIPPED_ALLELE_1, 
            :ALLELE_2 => :UNFLIPPED_ALLELE_2
        )
        flipped_variants_list = readlines(flip_file)
        flipped_variants = innerjoin(
            innerjoin(
                flipped_bim, 
                DataFrame(VARIANT_ID=flipped_variants_list),
                on=:VARIANT_ID
            ),
            unflipped_bim, 
            on=:VARIANT_ID
        )
        @test nrow(flipped_variants) == length(flipped_variants_list)
        complement = Dict("A" => "T", "T" => "A", "C" => "G", "G" => "C")
        for row in eachrow(flipped_variants)
            @test row.ALLELE_1 == complement[row.UNFLIPPED_ALLELE_1]
            @test row.ALLELE_2 == complement[row.UNFLIPPED_ALLELE_2]
        end
    end
    ## New fam files
    release_2024_now_fam_before_extract = SequentialGWAS.read_fam(joinpath(qced_dir, "mock.release_2024_now.qced.fam"))
    release_2024_now_fam_after_extract = SequentialGWAS.read_fam(joinpath(flipped_shared_dir, "release-2024-now.flipped.shared.fam"))
    release_2024_now_samples_to_drop = CSV.read(joinpath(kgp_qc_files_dir, "mock.release_2024_now.qced.samples_to_drop.txt"), DataFrame, header=false)
    dropped_samples_from_fam = setdiff(release_2024_now_fam_before_extract.IID, release_2024_now_fam_after_extract.IID)
    @test dropped_samples_from_fam == release_2024_now_samples_to_drop[!, 2] == ["odap3001"]

    release_r8_fam_before_extract = SequentialGWAS.read_fam(joinpath(qced_dir, "mock.release_r8.liftedOver.qced.fam"))
    release_r8_fam_after_extract = SequentialGWAS.read_fam(joinpath(flipped_shared_dir, "release-r8.flipped.shared.fam"))
    release_r8_samples_to_drop = CSV.read(joinpath(kgp_qc_files_dir, "mock.release_r8.liftedOver.qced.samples_to_drop.txt"), DataFrame, header=false)
    dropped_samples_from_fam = setdiff(release_r8_fam_before_extract.IID, release_r8_fam_after_extract.IID)
    @test sort(dropped_samples_from_fam) == sort(release_r8_samples_to_drop[!, 2]) == ["odap2002", "odap3001"]

    # Check WGS
    wgs_dir = joinpath(RESULTS_DIR, "wgs", "genotyped")
    gatk_shared_variants = CSV.read(
        joinpath(kgp_qc_files_dir, "variants_intersection.bed"), 
        DataFrame, 
        header=[:CHR, :BP_START, :BP_END]
    )
    @test length(shared_variants) == nrow(gatk_shared_variants)
    @test all(gatk_shared_variants.BP_START == gatk_shared_variants.BP_END .- 1)
    @test all(gatk_shared_variants.BP_END == [parse(Int, x[2]) for x in split.(shared_variants, ":")])

    odap_bim_files = joinpath.(wgs_dir, string.("mock.odap", 3001:3010, ".shared.bim"))
    for bim_file in odap_bim_files
        bim = SequentialGWAS.read_bim(bim_file)
        # Check chr:pos:ref:alt format
        @test all(length.(split.(bim.VARIANT_ID, ":")) .== 4)
        # Check exactly required variants were genotyped
        @test length(setdiff(shared_variants, bim.VARIANT_ID)) == 0
        @test length(setdiff(bim.VARIANT_ID, shared_variants)) == 0
    end
    
    # Check merged genotypes
    merge_dir = joinpath(RESULTS_DIR, "merged")
    merged_bim = SequentialGWAS.read_bim(joinpath(merge_dir, "merged", "genotypes.merged.bim"))
    @test sort(merged_bim.VARIANT_ID) == sort(shared_variants)
    merged_fam = SequentialGWAS.read_fam(joinpath(merge_dir, "merged", "genotypes.merged.fam"))

    ## Check plink resolves inconsistent alleles when Merging
    ## Even though according to the following link this is the case, we make sure of this: https://wanggroup.org/compbio_tutorial/allele_qc.html
    ## for instance this variant has REF == G on the r8_release/2021_2023_release and REF == A on the 2024_release
    non_consistent_alleles_snp = "chr1:18100537:G:A"
    release_r8_freqs = CSV.read(joinpath(flipped_shared_dir, "release-r8.flipped.shared.acount"), DataFrame)
    release_2021_2023_freqs = CSV.read(joinpath(flipped_shared_dir, "release-2021-2023.flipped.shared.acount"), DataFrame)
    release_2024_now_freqs = CSV.read(joinpath(flipped_shared_dir, "release-2024-now.flipped.shared.acount"), DataFrame)
    wgs_freqs = readdlm(joinpath(RESULTS_DIR, "wgs", "merged", "wgs.merged.frq.counts"))
    wgs_freqs = DataFrame(wgs_freqs[2:end, :], Symbol.(wgs_freqs[1, :]))
    merged_freqs = readdlm(joinpath(merge_dir, "merged", "genotypes.merged.frq.counts"))
    merged_freqs = DataFrame(merged_freqs[2:end, :], Symbol.(merged_freqs[1, :]))

    release_r8_stats = only(filter(:ID => ==(non_consistent_alleles_snp), release_r8_freqs))
    release_2021_2023_stats = only(filter(:ID => ==(non_consistent_alleles_snp), release_2021_2023_freqs))
    release_2024_now_stats = only(filter(:ID => ==(non_consistent_alleles_snp), release_2024_now_freqs))
    wgs_stats = only(filter(:SNP => ==(non_consistent_alleles_snp), wgs_freqs))
    merged_stats = only(filter(:SNP => ==(non_consistent_alleles_snp), merged_freqs))

    ## Let's count the sum of A before and after merging
    @test release_r8_stats.REF == "G" == release_2021_2023_stats.REF !== release_2024_now_stats.REF == "A" == wgs_stats.A1
    @test release_r8_stats.ALT_CTS + release_2021_2023_stats.ALT_CTS + (release_2024_now_stats.OBS_CT - release_2024_now_stats.ALT_CTS) + wgs_stats.C1 == merged_stats.C1


    # Check QC of merged genotypes
    ## Check filtered samples
    qced_merged_bim = SequentialGWAS.read_bim(joinpath(merge_dir, "qced", "genotypes.merged.qced.bim"))
    @test issubset(qced_merged_bim.VARIANT_ID, shared_variants)
    qced_merged_fam = SequentialGWAS.read_fam(joinpath(merge_dir, "qced", "genotypes.merged.qced.fam"))
    unrelated_king = CSV.read(
        joinpath(RESULTS_DIR, "merged/king_relatedness/kingunrelated.txt"),
        DataFrame,
        header=[:INDEX, :IID]
    )
    @test nrow(qced_merged_fam) < nrow(merged_fam)
    @test issubset(qced_merged_fam.IID, unrelated_king.IID)

    ## Check Ancestry
    ancestry_dir = joinpath(RESULTS_DIR, "ancestry")
    kgp_merged_fam = SequentialGWAS.read_fam(joinpath(ancestry_dir, "merged","genotypes_and_kgp.merged.fam"))
    kgp_ancestries = CSV.read(
        joinpath(TESTDIR, "assets", "kgp", "20130606_g1k_3202_samples_ped_population.txt"), 
        DataFrame,
        select=[:SampleID, :Superpopulation]
    )
    leftjoin!(kgp_merged_fam, kgp_ancestries, on=:IID => :SampleID)
    pop_counts = combine(groupby(kgp_merged_fam, :Superpopulation), nrow)
    @test nrow(pop_counts) == 6
    Q = readdlm(joinpath(ancestry_dir, "ancestry", "genotypes_and_kgp.merged.ldpruned.5.Q"))
    kgp_merged_fam.POP_IDX = [i.I[2] for i in argmax(Q, dims=2)[:, 1]]
    kgp_fam = filter(x -> x.Superpopulation !== missing, kgp_merged_fam)
    # Check the prediction is the same individuals with known superpopulation
    @test length(groupby(kgp_fam, [:Superpopulation, :POP_IDX])) == 5
    # Check predictions are consistent
    ancestry_predictions = CSV.read(
        joinpath(ancestry_dir, "ancestry", "genotypes_and_kgp.merged.ldpruned.ancestry.csv"), 
        DataFrame, 
        select=[:IID, :Superpopulation]
    )
    genomicc_with_pred = innerjoin(
        kgp_merged_fam, 
        select(ancestry_predictions, :IID, :Superpopulation => :SuperpopulationPred),
        on=:IID
    )
    @test length(groupby(genomicc_with_pred, [:SuperpopulationPred, :POP_IDX])) > 0
    
    # Check PCA-QC
    pca_dir = joinpath(merge_dir, "pca_qced")
    @test isfile(joinpath(pca_dir, "genotypes.merged.qced.ldpruned.after_pca_qc.1vs2.png"))
    @test isfile(joinpath(pca_dir, "genotypes.merged.qced.ldpruned.after_pca_qc.all.png"))
    @test isfile(joinpath(pca_dir, "genotypes.merged.qced.ldpruned.before_pca_qc.loadings.png"))
    @test isfile(joinpath(pca_dir, "genotypes.merged.qced.ldpruned.after_pca_qc.loadings.png"))

    # Check report and final dataset
    @test isfile(joinpath(RESULTS_DIR, "report.md"))
    @test isfile(joinpath(RESULTS_DIR, "genotypes.arrays_wgs.aggregated.bed"))
    @test isfile(joinpath(RESULTS_DIR, "genotypes.arrays_wgs.aggregated.bim"))
    @test isfile(joinpath(RESULTS_DIR, "genotypes.arrays_wgs.aggregated.fam"))
    final_bim = SequentialGWAS.read_bim(joinpath(RESULTS_DIR, "genotypes.arrays_wgs.aggregated.bim"))
    variants_removed_by_pca = readlines(joinpath(pca_dir, "variants_to_exclude.txt"))
    @test intersect(final_bim.VARIANT_ID, variants_removed_by_pca) == []

end

end
true