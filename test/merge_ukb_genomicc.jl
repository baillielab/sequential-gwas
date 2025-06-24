module TestUKBGenOMICCMerge

using Test
using SequentialGWAS
using BGEN
using CSV
using DataFrames

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

        results_dirs = readdir("ukb_genomicc_merge_results/merge_ukb_and_genomicc/", join=true)
        results_dir = results_dirs[argmax(mtime(d) for d in results_dirs)]

        # Test filtering chromosome files (3 filesets)
        all_chr_bim = DataFrame()
        for shard in [0, 1, 2]
            execution_dir = joinpath(results_dir, "call-filter_ukb_chr", "shard-$shard", "execution")
            results_files = readdir(execution_dir, join=true)
            # Read bim file
            bim_file = only(filter(x -> endswith(x, ".bim"), results_files))
            bim = SequentialGWAS.read_bim(bim_file)
            append!(all_chr_bim, DataFrame(bim))
            chr = replace(first(bim.CHR_CODE), "chr" => "")
            # Read fam file
            fam_file = only(filter(x -> endswith(x, ".fam"), results_files))
            fam = SequentialGWAS.read_fam(fam_file)
            # Check variants are subset of source BGEN file
            src_bgen_file = joinpath(TESTDIR, "assets", "ukb", "imputed", "ukb21007_c$(chr)_b0_v1.bgen")
            b = BGEN.Bgen(src_bgen_file)
            src_positions = [pos(v) for v in iterator(b)]
            @test issubset(bim.BP_COORD, src_positions)
            # Check critical care individuals have been dropped
            samples_in_critical_care = CSV.read(
                joinpath(TESTDIR, "assets", "ukb", "critical_table.csv"), 
                DataFrame
            ).eid
            samples_in_filterd_chr = fam.IID
            @test isempty(intersect(fam.IID, samples_in_critical_care))
        end

        # Test merging chromosome files
        merged_ukb_chr_dir = joinpath(results_dir, "call-merge_ukb_chrs", "execution")
        merged_chr_bim = SequentialGWAS.read_bim(joinpath(merged_ukb_chr_dir, "ukb_all_chr.bim"))
        @test sort(merged_chr_bim) == sort(all_chr_bim)
        merged_chr_fam = SequentialGWAS.read_fam(joinpath(merged_ukb_chr_dir, "ukb_all_chr.fam"))

        # Test variant Ids alignement and filterting of unrelated individuals
        kgp_qc_dir = joinpath(results_dir, "call-align_ukb_variant_ids_with_kgp_and_keep_unrelated", "execution")
        ukb_unrelated_bim = SequentialGWAS.read_bim(joinpath(kgp_qc_dir, "ukb_unrelated.bim"))
        ## One variant is dropped because it is not in the KGP dataset
        @test 9694126 ∉ ukb_unrelated_bim.BP_COORD
        ## Samples, none is removed here.
        ukb_unrelated_fam = SequentialGWAS.read_fam(joinpath(kgp_qc_dir, "ukb_unrelated.fam"))
        @test length(ukb_unrelated_fam.IID) <= length(merged_chr_fam.IID)
        
        # Test merging UKB and KGP datasets
        ukb_kgp_merged_dir = joinpath(results_dir, "call-merge_ukb_kgp", "execution")
        merged_ukb_kgp_bim = SequentialGWAS.read_bim(joinpath(ukb_kgp_merged_dir, "ukb_kgp.merged.bim"))
        @test length(merged_ukb_kgp_bim.VARIANT_ID) <= length(ukb_unrelated_bim.VARIANT_ID)
        merged_ukb_kgp_fam = SequentialGWAS.read_fam(joinpath(ukb_kgp_merged_dir, "ukb_kgp.merged.fam"))
        @test length(merged_ukb_kgp_fam.IID) > length(ukb_unrelated_fam.IID)

        # Test LD pruning
        ld_pruning_dir = joinpath(results_dir, "call-ld_prune_ukb_kgp", "execution")
        ld_pruned_bim = SequentialGWAS.read_bim(joinpath(ld_pruning_dir, "ukb_kgp.merged.qc.ld_pruned.bim"))
        @test length(ld_pruned_bim.VARIANT_ID) <= length(merged_ukb_kgp_bim.VARIANT_ID)

        # Test ancestry estimation
        ancestry_file = joinpath(results_dir, "call-estimate_ukb_ancestry_from_kgp", "execution", "ukb.ancestry_estimate.csv")
        ancestry_df = CSV.read(ancestry_file, DataFrame)
        @test length(ancestry_df.IID) == length(ukb_unrelated_fam.IID)
        @test names(ancestry_df) == ["FID", "IID", "Superpopulation", "AFR", "SAS", "EAS", "AMR", "EUR"]

        # Test merging with GenOMICC
        ukb_genomicc_merged_dir = joinpath(results_dir, "call-merge_ukb_genomicc", "execution")
        ukb_genomicc_merged_bim = SequentialGWAS.read_bim(joinpath(ukb_genomicc_merged_dir, "ukb_genomicc.merged.bim"))
        @test length(ukb_genomicc_merged_bim.VARIANT_ID) <= length(ukb_unrelated_bim.VARIANT_ID)
        ukb_genomicc_merged_fam = SequentialGWAS.read_fam(joinpath(ukb_genomicc_merged_dir, "ukb_genomicc.merged.fam"))
        @test length(ukb_genomicc_merged_fam.IID) > length(ukb_unrelated_fam.IID)
    end
end
end

true