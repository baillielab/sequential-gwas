module TestUKBGenOMICCMerge

using Test
using GenomiccWorkflows
using BGEN
using CSV
using DataFrames

PKGDIR = pkgdir(GenomiccWorkflows)
TESTDIR = joinpath(PKGDIR, "test")

@testset "Test misc functions" begin
    # Test format_chromosome! function
    bim = DataFrame(
        CHR_CODE=["chr1", "chr2", "chr3"],
        BP_COORD=[1000, 2000, 3000],
        VARIANT_ID=["rs1", "rs2", "rs3"]
    )
    GenomiccWorkflows.format_chromosome!(bim)
    @test all(bim.CHR_CODE .== ["1", "2", "3"])
    
    bim = DataFrame(
        CHR_CODE=[1, 2, 3],
        BP_COORD=[1000, 2000, 3000],
        VARIANT_ID=["rsX1", "rsY1", "rsM1"],
        ALLELE_1=["A", "G", "C"],
        ALLELE_2=["C", "T", "G"]
    )
    GenomiccWorkflows.format_chromosome!(bim)
    @test all(bim.CHR_CODE .== ["1", "2", "3"])

    # Test update_variant_ids_with_map! function
    variant_ids_map = Dict(
        ("1", 1000) => ("rs1_kgp", Set(["A", "C"])), 
        ("2", 2000) => ("rs2_kgp", Set(["G", "C"]))
    )
    unmapped_ids = GenomiccWorkflows.update_variant_ids_with_map!(bim, variant_ids_map)
    @test all(bim.VARIANT_ID .== ["rs1_kgp", "rsY1", "rsM1"])
    @test unmapped_ids == Set(["rsM1", "rsY1"])
end

@testset "Test merge_ukb_genomicc_covariates" begin
    tmpdir = mktempdir()
    output_file = joinpath(tmpdir, "ukb_genomicc.covariates.csv")
    genomicc_covariates_file = joinpath("test", "assets", "genomicc", "mock.covariates.csv")
    genomicc_inferred_covariates_file = joinpath("test", "assets", "genomicc", "inferred_covariates.csv")
    ukb_covariates_file = joinpath("test", "assets", "ukb", "covariates_table.csv")
    ukb_inferred_covariates_file = joinpath(tmpdir, "inferred_covariates.ukb.csv")
    file_with_eids_to_exclude = joinpath(TESTDIR, "assets", "ukb", "critical_table.csv")
    # Make UKB inferred covariates (ancestry estimates and downsample to mimic loss of samples)
    ukb_covariates = CSV.read(ukb_covariates_file, DataFrame)
    n = nrow(ukb_covariates)
    ukb_inferred_covariates = DataFrame(
        FID = ukb_covariates.eid,
        IID = ukb_covariates.eid,
        Superpopulation = rand(["AFR", "SAS", "EAS", "AMR", "EUR"], n),
        AFR = rand(n),
        SAS = rand(n),
        EAS = rand(n),
        AMR = rand(n),
        EUR = rand(n)
    )
    ukb_inferred_covariates = ukb_inferred_covariates[1:15, :]
    CSV.write(ukb_inferred_covariates_file, ukb_inferred_covariates)
    # Merge covariates
    copy!(ARGS, [
        "merge-ukb-genomicc-covariates",
        genomicc_covariates_file,
        genomicc_inferred_covariates_file,
        ukb_covariates_file,
        ukb_inferred_covariates_file,
        file_with_eids_to_exclude,
        "--output-file", output_file
    ])
    julia_main()
    merged_covariates = CSV.read(output_file, DataFrame)
    @test nrow(merged_covariates) == 14 + 12000 # ukb19 is dropped because in critical_table.csv
    @test names(merged_covariates) == ["FID", "IID", "AGE", "SEX", "ANCESTRY_ESTIMATE", "AFR", "SAS", "EAS", "AMR", "EUR", "COHORT"]
    @test Set(merged_covariates.COHORT) == Set(["GENOMICC", "UKB"])
    @test eltype(merged_covariates.AGE) == Int
    @test Set(merged_covariates.SEX) == Set([0, 1, missing])
    @test Set(merged_covariates.ANCESTRY_ESTIMATE) == Set(["AFR", "SAS", "EAS", "AMR", "EUR"])
    for col in [:AFR, :SAS, :EAS, :AMR, :EUR]
        @test eltype(merged_covariates[!, col]) == Float64
    end
    @test merged_covariates.FID == merged_covariates.IID
end


@testset "Test make_ukb_bgen_qc_and_r2_filter_files" begin
    tmpdir = mktempdir()
    prefix = joinpath(TESTDIR, "assets", "ukb", "imputed", "ukb21007_c1_b0_v1")
    output_file = joinpath(tmpdir, "extract_list.txt")
    pvar = CSV.read(string(prefix, ".pvar"), DataFrame; delim='\t')
    copy!(ARGS, [
        "make-ukb-bgen-qc-and-r2-filter-files",
        prefix,
        "--threshold", "0.9",
        "--output", output_file
    ])
    julia_main()
    # First variant with R2 < 0.9 is dropped
    extract_list = readlines(output_file)
    @test extract_list == [
        "chr1:14012312:T:C",
        "chr1:18100537:G:A",
        "chr1:22542609:T:C",
        "chr1:40310265:G:A",
        "chr1:92682820:C:T",
        "chr1:111622622:C:A",
        "chr1:183905563:G:A",
        "chr1:231799576:C:T"
    ]
    # pvar file has updated ID
    new_pvar = CSV.read(string(prefix, ".pvar"), DataFrame; delim='\t')
    @test new_pvar.ID == [
        "chr1:9694126:C:T",
        "chr1:14012312:T:C",
        "chr1:18100537:G:A",
        "chr1:22542609:T:C",
        "chr1:40310265:G:A",
        "chr1:92682820:C:T",
        "chr1:111622622:C:A",
        "chr1:183905563:G:A",
        "chr1:231799576:C:T"
    ]
    # Write back the original pvar file
    CSV.write(
        string(prefix, ".pvar"),
        pvar,
        delim='\t',
        writeheader=true
    )
end

# End to End Workflow run

dorun = isinteractive() || (haskey(ENV, "CI_CONTAINER") && ENV["CI_CONTAINER"] == "docker")

if dorun
    cmd_args = haskey(ENV, "CROMWELL_PATH") ?
        ["-jar", ENV["CROMWELL_PATH"]] :
        ["-Dconfig.file=conf/cromwell.mac.conf", "-jar", "/Users/olabayle/cromwell/cromwell-90.jar"]

    # Download the reference genome if not already present
    reference_genome_path = joinpath(PKGDIR, "assets", "rap", "Homo_sapiens_assembly38.fasta")
    if !isfile(reference_genome_path)
        run(`wget -O $reference_genome_path https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta`)
    end

    @testset "Test Array Genotypes Merging" begin
        cmd = Cmd([
            "java", cmd_args...,
            "run", "rap_workflows/ukb_merge/workflow.wdl",
            "--inputs", joinpath(TESTDIR, "assets", "ukb_merge.json"),
            "--options", joinpath(TESTDIR, "assets", "ukb_merge_options.json")
        ])
        rc = run(cmd)
        @test rc.exitcode == 0

        results_dirs = readdir("ukb_genomicc_merge_results/merge_ukb_and_genomicc/", join=true)
        results_dir = results_dirs[argmax(mtime(d) for d in results_dirs)]

        # Test R2 filter and critical samples removed
        samples_in_critical_care = CSV.read(
                joinpath(TESTDIR, "assets", "ukb", "critical_table.csv"), 
                DataFrame
            ).eid
        first_filter_dir = joinpath(results_dir, "call-filter_ukb_chr_with_r2_and_critical_samples")
        for shard in [0, 1, 2]
            execution_dir = joinpath(first_filter_dir, "shard-$shard", "execution")
            results_files = readdir(execution_dir, join=true)
            # Check critical samples are removed
            samples = CSV.read(
                only(filter(x -> endswith(x, ".psam"), results_files)),
                DataFrame
            )
            @test isempty(intersect(samples.IID, samples_in_critical_care))
            # Check R2 filter is applied
            pvar = CSV.read(
                only(filter(x -> endswith(x, ".pvar"), results_files)),
                DataFrame
            )
            if first(pvar[!, "#CHROM"]) == "chr1"
                @test "chr1:9694126:C:T" ∉ pvar.ID
            end
        end

        # Test extract genomicc variants
        all_chr_bim = DataFrame()
        for shard in [0, 1, 2]
            execution_dir = joinpath(results_dir, "call-extract_genomicc_variants", "shard-$shard", "execution")
            results_files = readdir(execution_dir, join=true)
            # Read bim file
            bim_file = only(filter(x -> endswith(x, ".bim"), results_files))
            bim = GenomiccWorkflows.read_bim(bim_file)
            append!(all_chr_bim, DataFrame(bim))
        end

        # Test merging chromosome files
        merged_ukb_chr_dir = joinpath(results_dir, "call-merge_ukb_chrs", "execution")
        merged_chr_bim = GenomiccWorkflows.read_bim(joinpath(merged_ukb_chr_dir, "ukb_all_chr.bim"))
        @test sort(merged_chr_bim) == sort(all_chr_bim)
        merged_chr_fam = GenomiccWorkflows.read_fam(joinpath(merged_ukb_chr_dir, "ukb_all_chr.fam"))

        # Test variant Ids alignement and filterting of unrelated individuals
        kgp_qc_dir = joinpath(results_dir, "call-align_ukb_variants_with_kgp_and_keep_unrelated", "execution")
        ukb_unrelated_bim = GenomiccWorkflows.read_bim(joinpath(kgp_qc_dir, "ukb_unrelated.bim"))
        ## One variant is dropped because it is not in the KGP dataset
        @test 9694126 ∉ ukb_unrelated_bim.BP_COORD
        ## Samples, none is removed here.
        ukb_unrelated_fam = GenomiccWorkflows.read_fam(joinpath(kgp_qc_dir, "ukb_unrelated.fam"))
        @test length(ukb_unrelated_fam.IID) <= length(merged_chr_fam.IID)
        
        # Test merging UKB and KGP datasets
        ukb_kgp_merged_dir = joinpath(results_dir, "call-merge_ukb_kgp", "execution")
        merged_ukb_kgp_bim = GenomiccWorkflows.read_bim(joinpath(ukb_kgp_merged_dir, "ukb_kgp.merged.bim"))
        @test length(merged_ukb_kgp_bim.VARIANT_ID) <= length(ukb_unrelated_bim.VARIANT_ID)
        merged_ukb_kgp_fam = GenomiccWorkflows.read_fam(joinpath(ukb_kgp_merged_dir, "ukb_kgp.merged.fam"))
        @test length(merged_ukb_kgp_fam.IID) > length(ukb_unrelated_fam.IID)

        # Test LD pruning
        ld_pruning_dir = joinpath(results_dir, "call-ld_prune_ukb_kgp", "execution")
        ld_pruned_bim = GenomiccWorkflows.read_bim(joinpath(ld_pruning_dir, "ukb_kgp.merged.qc.ld_pruned.bim"))
        @test length(ld_pruned_bim.VARIANT_ID) <= length(merged_ukb_kgp_bim.VARIANT_ID)

        # Test ancestry estimation
        ancestry_file = joinpath(results_dir, "call-estimate_ukb_ancestry_from_kgp", "execution", "ukb.ancestry_estimate.csv")
        ancestry_df = CSV.read(ancestry_file, DataFrame)
        @test length(ancestry_df.IID) == length(ukb_unrelated_fam.IID)
        @test names(ancestry_df) == ["FID", "IID", "Superpopulation", "AFR", "SAS", "EAS", "AMR", "EUR"]

        # Test merging with GenOMICC
        ukb_genomicc_merged_dir = joinpath(results_dir, "call-merge_ukb_genomicc", "execution")
        ukb_genomicc_merged_bim = GenomiccWorkflows.read_bim(joinpath(ukb_genomicc_merged_dir, "ukb_genomicc.merged.bim"))
        @test length(ukb_genomicc_merged_bim.VARIANT_ID) <= length(ukb_unrelated_bim.VARIANT_ID)
        ukb_genomicc_merged_fam = GenomiccWorkflows.read_fam(joinpath(ukb_genomicc_merged_dir, "ukb_genomicc.merged.fam"))
        @test length(ukb_genomicc_merged_fam.IID) > length(ukb_unrelated_fam.IID)

        # Test merging covariates
        ukb_genomicc_covariates_file = joinpath(results_dir, "call-merge_genomicc_ukb_covariates", "execution", "ukb_genomicc.covariates.csv")
        ukb_genomicc_covariates = CSV.read(ukb_genomicc_covariates_file, DataFrame)
        @test names(ukb_genomicc_covariates) == ["FID", "IID", "AGE", "SEX", "ANCESTRY_ESTIMATE", "AFR", "SAS", "EAS", "AMR", "EUR", "COHORT"]
        @test length(filter(startswith("ukb"), ukb_genomicc_covariates.IID)) > 0
        @test length(filter(startswith("odap"), ukb_genomicc_covariates.IID)) > 0

        # Test report generation
        @test isfile(joinpath(results_dir, "call-make_report", "execution", "report.md"))
    end
end

end

true