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

    # Test add_primary_diagnosis!
    odap_covariates = DataFrame(PRIM_DIAGNOSIS_ODAP=[
        "covid-19", "isaric4c covid-19", "mild covid-19", "react covid-19",
        "influenza virus", "pneumonia with radiographic changes at presentation to critical care",
        "pancreatitis of any aetiology", "rsv (respiratory syncytial virus) infection",
        "soft tissue infections causing systemic sepsis", "ecls", "reaction to vaccination"
    ])
    GenomiccWorkflows.add_primary_diagnosis!(odap_covariates)
    @test odap_covariates.PRIMARY_DIAGNOSIS == [
        "COVID_19", "COVID_19", "COVID_19", "COVID_19",
        "INFLUENZA", "PNEUMONIA",
        "PANCREATITIS", "RSV",
        "SOFT_TISSUE_INFECTION", "ECLS", "REACTION_TO_VACCINATION"
    ]
    odap_covariates = DataFrame(PRIM_DIAGNOSIS_ODAP=["unknown diagnosis"])
    @test_throws KeyError GenomiccWorkflows.add_primary_diagnosis!(odap_covariates)

    # Test is_severe_covid_19
    covariates = DataFrame(
        PRIMARY_DIAGNOSIS=["COVID_19", "flu", "COVID_19", "COVID_19", "COVID_19", "COVID_19", "COVID_19", "COVID_19"],
        COHORT=["GENOMICC_SEVERE", "GENOMICC_SEVERE", "GEN_INT_PAKISTAN", "GENOMICC_MILD", "GENOMICC_REACT", "ISARIC4C", "ISARIC4C", "ISARIC4C"],
        ISARIC_MAX_SEVERITY_SCORE=[missing, missing, missing, missing, missing, "NA", "3", "4"],
        EXPECTED_SEVERELY_ILL = [1, 1, 1, 0, 0, missing, 0, 1],
        EXPECTED_SEVERE_COVID_19=[1, missing, 1, 0, 0, missing, 0, 1]
    )
    GenomiccWorkflows.add_is_severely_ill_col!(covariates)
    @test all(covariates.IS_SEVERELY_ILL .=== covariates.EXPECTED_SEVERELY_ILL)
    covariates.SEVERE_COVID_19 = map(x -> GenomiccWorkflows.is_severe_infection(x, "COVID_19"), eachrow(covariates))
    @test all(covariates.SEVERE_COVID_19 .=== covariates.EXPECTED_SEVERE_COVID_19)
    unknown_cohort_covariates = DataFrame(
        PRIMARY_DIAGNOSIS=["COVID_19"],
        COHORT=["unknown cohort"],
        ISARIC_MAX_SEVERITY_SCORE=[missing]
    )
    @test_throws ArgumentError("Unknown cohort: unknown cohort") GenomiccWorkflows.add_is_severely_ill_col!(unknown_cohort_covariates)
    
    # Test is_severe_infection_influenza
    covariates = DataFrame(
        PRIMARY_DIAGNOSIS=["COVID_19", "INFLUENZA", "INFLUENZA"],
        ISARIC_MAX_SEVERITY_SCORE=["5", missing, missing],
        COHORT=["GENOMICC_SEVERE", "GENOMICC_SEVERE", "GEN_INT_PAKISTAN"],
        EXPECTED_IS_SEVERELY_ILL=[1, 1, 1],
        EXPECTED_SEVERE_INFLUENZA=[missing, 1, 1]
    )
    GenomiccWorkflows.add_is_severely_ill_col!(covariates)
    @test all(covariates.IS_SEVERELY_ILL .=== covariates.EXPECTED_IS_SEVERELY_ILL)
    covariates.SEVERE_INFLUENZA = map(x -> GenomiccWorkflows.is_severe_infection(x, "INFLUENZA"), eachrow(covariates))
    @test all(covariates.SEVERE_INFLUENZA .=== covariates.EXPECTED_SEVERE_INFLUENZA)
    # Test is_alive_at_assessment
    df = DataFrame(
        ALIVE_AT_60_DAYS = ["Yes", "no", "yes", "No", "NA", "NA", "NA", "NA", "NA"],
        ALIVE_AT_28_DAYS = ["NA", "NA", "NA", "NA", "yes", "Yes", "No", "no", "NA"]
    )
    GenomiccWorkflows.add_alive_at_assessment_col!(df)
    @test all(df.ALIVE_AT_ASSESSMENT .=== [1, 0, 1, 0, 1, 1, 0, 0, missing])
end

@testset "Test process_genomicc_covariates: without UKB" begin
    tmpdir = mktempdir()
    output_file = joinpath(tmpdir, "covariates.processed.csv")
    genomicc_covariates_file = joinpath(TESTDIR, "assets", "genomicc", "mock.covariates.csv")
    copy!(ARGS, [
        "process-genomicc-covariates",
        genomicc_covariates_file,
        "--output-file", output_file
    ])
    julia_main()
    processed_covariates = CSV.read(joinpath(tmpdir, "covariates.processed.csv"), DataFrame)
    @test Set(names(processed_covariates)) == Set([
        "FID",
        "IID",
        "COHORT",
        "PRIMARY_DIAGNOSIS",
        "AGE",
        "SEX",
        "SUPERPOPULATION",
        "AFR",
        "AMR",
        "EAS",
        "EUR",
        "SAS",
        "IS_SEVERELY_ILL",
        "ALIVE_AT_ASSESSMENT",
        "SEVERE_COVID_19",
        "SEVERE_PNEUMONIA",
        "SEVERE_PANCREATITIS",
        "SEVERE_INFLUENZA",
        "SEVERE_SOFT_TISSUE_INFECTION",
        "SEVERE_RSV",
        "SEVERE_ECLS",
        "SEVERE_REACTION_TO_VACCINATION"])
    @test countlines(output_file) == countlines(genomicc_covariates_file) # No samples are dropped
end

@testset "Test process_genomicc_covariates: with UKB" begin
    tmpdir = mktempdir()
    output_file = joinpath(tmpdir, "ukb_genomicc.covariates.csv")
    genomicc_covariates_file = joinpath(TESTDIR, "assets", "genomicc", "mock.covariates.csv")
    ukb_covariates_file = joinpath(TESTDIR, "assets", "ukb", "covariates_table.csv")
    ukb_inferred_covariates_file = joinpath(tmpdir, "inferred_covariates.ukb.csv")
    unrelated_ukb_individuals_file = joinpath(tmpdir, "kingunrelated.txt")
    # Make UKB inferred covariates (ancestry estimates and downsample to mimic loss of samples due to both relatedness and critical patients)
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
        "process-genomicc-covariates",
        genomicc_covariates_file,
        string("--ukb-covariates=", ukb_covariates_file),
        string("--ukb-inferred-covariates=", ukb_inferred_covariates_file),
        "--output-file", output_file
    ])
    julia_main()
    merged_covariates = CSV.read(output_file, DataFrame)
    ukb_individuals = merged_covariates.IID[merged_covariates.IID .|> x -> startswith(x, "ukb")]
    @test nrow(merged_covariates) == 15 + 12911 # ukb19 is dropped because in critical_table.csv
    @test Set(names(merged_covariates)) == Set([
        "FID",
        "IID",
        "COHORT",
        "PRIMARY_DIAGNOSIS",
        "AGE",
        "SEX",
        "SUPERPOPULATION",
        "AFR",
        "AMR",
        "EAS",
        "EUR",
        "SAS",
        "IS_SEVERELY_ILL",
        "ALIVE_AT_ASSESSMENT",
        "SEVERE_PNEUMONIA",
        "SEVERE_COVID_19",
        "SEVERE_PANCREATITIS",
        "SEVERE_INFLUENZA",
        "SEVERE_SOFT_TISSUE_INFECTION",
        "SEVERE_RSV",
        "SEVERE_ECLS",
        "SEVERE_REACTION_TO_VACCINATION"
    ])
    @test Set(merged_covariates.COHORT) == Set(["ISARIC4C", "GEN_INT_PAKISTAN", "GENOMICC_MILD", "GENOMICC_SEVERE", "GENOMICC_REACT", "UKBIOBANK"])
    @test eltype(merged_covariates.AGE) == Int # No missing value in AGE
    @test Set(merged_covariates.SEX) == Set([0, 1, missing]) # Individuals with missing Sex will be dropped in analyses
    @test Set(merged_covariates.SUPERPOPULATION) == Set(["AFR", "SAS", "EAS", "AMR", "EUR", "ADMIXED"])
    @test merged_covariates.FID == merged_covariates.IID
    @test Set(merged_covariates.ALIVE_AT_ASSESSMENT) == Set([0, 1, missing])
    @test Set(merged_covariates.PRIMARY_DIAGNOSIS) == Set([
        "COVID_19", "INFLUENZA", "PANCREATITIS", "PNEUMONIA",
        "RSV", "SOFT_TISSUE_INFECTION", "ECLS", "REACTION_TO_VACCINATION",
        missing
    ])
    @test all(merged_covariates[merged_covariates.COHORT .== "UKBIOBANK", :IS_SEVERELY_ILL] .== 0)
end

@testset "Test fill_chr_pvar_with_variant_id" begin
    tmpdir = mktempdir()
    # The file is modified inplace so we need to copy it to preserve the test dir
    pvar_file = joinpath(tmpdir, "temp.pvar")
    src_pvar_file = joinpath(TESTDIR, "assets", "ukb", "imputed", "ukb21007_c1_b0_v1.pvar")
    cp(src_pvar_file, pvar_file)
    variants_info_file = joinpath(TESTDIR, "assets", "ukb", "imputed", "ukb21007_c1_b0_v1.tsv")
    
    copy!(ARGS, [
        "fill-chr-pvar-with-variant-id",
        pvar_file, variants_info_file
    ])
    julia_main()
    # pvar file has updated ID
    src_pvar = CSV.read(src_pvar_file, DataFrame; delim='\t')
    @test all(src_pvar.ID .== ".")
    new_pvar = CSV.read(pvar_file, DataFrame; delim='\t')
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
end

@testset "Test make_ukb_individuals_list" begin
    tmpdir = mktempdir()
    covariates_file = joinpath(TESTDIR, "assets", "ukb", "covariates_table.csv")
    critical_table_file = joinpath(TESTDIR, "assets", "ukb", "critical_table.csv")
    original_eids = unique(CSV.read(covariates_file, DataFrame).eid)
    # Make UKB individuals list
    # This will filter out individuals with critical conditions
    output_file = joinpath(tmpdir, "ukb_eids_to_keep.txt")

    copy!(ARGS, [
        "make-ukb-individuals-list",
        covariates_file,
        critical_table_file,
        "--output", output_file,
    ])
    julia_main()

    eids_to_keep = readlines(output_file)
    @test length(eids_to_keep) == length(original_eids) - 5 # One individual (ukb19) is in critical_table.csv
    # Now further limit the number of individuals to 5
    copy!(ARGS, [
        "make-ukb-individuals-list",
        covariates_file,
        critical_table_file,
        "--output", output_file,
        "--max-samples", "5"
    ])
    julia_main()
    eids_to_keep = readlines(output_file)
    @test length(eids_to_keep) == 5
    @test "ukb19" ∉ eids_to_keep # ukb19 is still excluded
end

# End to End Workflow run

dorun = isinteractive() || (haskey(ENV, "CI_CONTAINER") && ENV["CI_CONTAINER"] == "docker")

if dorun
    # Download the reference genome if not already present
    reference_genome_path = joinpath(PKGDIR, "assets", "rap", "Homo_sapiens_assembly38.fasta")
    if !isfile(reference_genome_path)
        run(`wget -O $reference_genome_path https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta`)
    end

    @testset "Test Array Genotypes Merging" begin
        configfile = Sys.isapple() ? joinpath("conf", "cromwell.mac.conf") : joinpath("conf", "cromwell.local.conf")
        cmd = Cmd([
            "java", "-Dconfig.file=$configfile",
            "-jar", ENV["CROMWELL_PATH"],
            "run", joinpath(PKGDIR, "rap_workflows", "ukb_merge", "workflow.wdl"),
            "--inputs", joinpath(TESTDIR, "assets", "ukb_merge.json"),
            "--options", joinpath(TESTDIR, "assets", "ukb_merge_options.json")
        ])
        cd(PKGDIR) do
            rc = run(cmd)
            @test rc.exitcode == 0
        end

        results_dirs = readdir(joinpath(PKGDIR, "ukb_genomicc_merge_results", "merge_ukb_and_genomicc"), join=true)
        results_dir = results_dirs[argmax(mtime(d) for d in results_dirs)]

        # Test get_ukb_individuals
        ukb_individuals = CSV.read(
            joinpath(results_dir, "call-get_ukb_individuals", "execution", "ukb_eids_to_keep.txt"), 
            DataFrame;
            header=["FID", "IID"]
        )
        @test length(ukb_individuals.IID) == 150
        @test "ukb19" ∉ ukb_individuals.IID # ukb19 is excluded from the list

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
            @test length(samples.IID) <= 150 # missing rate may result in less than 15 individuals
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
        nindiv_per_files = []
        for shard in [0, 1, 2]
            execution_dir = joinpath(results_dir, "call-extract_genomicc_variants", "shard-$shard", "execution")
            results_files = readdir(execution_dir, join=true)
            # Read bim file
            bim_file = only(filter(x -> endswith(x, ".bim"), results_files))
            bim = GenomiccWorkflows.read_bim(bim_file)
            append!(all_chr_bim, DataFrame(bim))
            # Read fam file
            fam_file = only(filter(x -> endswith(x, ".fam"), results_files))
            fam = GenomiccWorkflows.read_fam(fam_file)
            push!(nindiv_per_files, nrow(fam))
        end
        ## Each imputed file dros some individuals, here we make sure only the intersection is kept for downstream analysis
        @test all(nindiv_per_files .== 145)

        # Test merging chromosome files
        merged_ukb_chr_dir = joinpath(results_dir, "call-merge_ukb_chrs", "execution")
        merged_chr_bim = GenomiccWorkflows.read_bim(joinpath(merged_ukb_chr_dir, "ukb_all_chr.bim"))
        @test sort(merged_chr_bim) == sort(all_chr_bim)
        merged_chr_fam = GenomiccWorkflows.read_fam(joinpath(merged_ukb_chr_dir, "ukb_all_chr.fam"))
        @test nrow(merged_chr_fam) == 145 # no more filtering of individuals

        # Test variant Ids alignement and filterting of unrelated individuals
        kgp_qc_dir = joinpath(results_dir, "call-align_ukb_variants_with_kgp_and_keep_unrelated", "execution")
        ukb_unrelated_bim = GenomiccWorkflows.read_bim(joinpath(kgp_qc_dir, "ukb_unrelated.bim"))
        unrelated_ukb_indivifuals = CSV.read(
            joinpath(kgp_qc_dir, "kingunrelated.txt"), 
            DataFrame; 
            delim="\t", 
            header=["FID", "IID"]
        )
        @test nrow(unrelated_ukb_indivifuals) <= 145

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
        @test names(ancestry_df) == ["FID", "IID", "Superpopulation", "AFR", "AMR", "EAS",  "EUR", "SAS"]

        # Test merging with GenOMICC
        ukb_genomicc_merged_dir = joinpath(results_dir, "call-merge_ukb_genomicc", "execution")
        ukb_genomicc_merged_bim = GenomiccWorkflows.read_bim(joinpath(ukb_genomicc_merged_dir, "ukb_genomicc.merged.bim"))
        @test length(ukb_genomicc_merged_bim.VARIANT_ID) <= length(ukb_unrelated_bim.VARIANT_ID)
        ukb_genomicc_merged_fam = GenomiccWorkflows.read_fam(joinpath(ukb_genomicc_merged_dir, "ukb_genomicc.merged.fam"))
        @test length(ukb_genomicc_merged_fam.IID) > length(ukb_unrelated_fam.IID)

        # Test merging covariates
        ukb_genomicc_covariates_file = joinpath(results_dir, "call-merge_genomicc_ukb_covariates", "execution", "ukb_genomicc.covariates.csv")
        ukb_genomicc_covariates = CSV.read(ukb_genomicc_covariates_file, DataFrame)
        @test Set(names(ukb_genomicc_covariates)) == Set([
            "FID", "IID", "COHORT", "AGE", "SEX", "SUPERPOPULATION",
            "AFR", "AMR", "EAS", "EUR", "SAS", "ALIVE_AT_ASSESSMENT",
            "SEVERE_PNEUMONIA", "SEVERE_COVID_19", "SEVERE_PANCREATITIS", 
            "SEVERE_INFLUENZA", "SEVERE_SOFT_TISSUE_INFECTION", "SEVERE_RSV", 
            "SEVERE_ECLS", "SEVERE_REACTION_TO_VACCINATION", "PRIMARY_DIAGNOSIS",
            "IS_SEVERELY_ILL"
        ])
        @test length(filter(startswith("ukb"), ukb_genomicc_covariates.IID)) > 0
        @test length(filter(startswith("odap"), ukb_genomicc_covariates.IID)) > 0

        # Imputed genotypes Merged
        ukb_genomicc_imputed_dir = joinpath(results_dir, "call-merge_genomicc_ukb_bcfs_and_convert_to_pgen")
        imputed_genotypes_individuals = []
        for shard in [0, 1, 2]
            subdir = joinpath(ukb_genomicc_imputed_dir, "shard-$shard", "execution")
            files = readdir(subdir)
            psam_file = files[findfirst(endswith(".psam"), files)]
            psam = CSV.read(joinpath(subdir, psam_file), DataFrame; delim='\t')
            rename!(psam, "#FID" => "FID")
            push!(imputed_genotypes_individuals, sort(psam[!, ["FID", "IID"]]))
        end

        # Test imputed genotypes, genotypes and genotypes have the same individuals
        for imputed_chr in imputed_genotypes_individuals
            @test imputed_chr == sort(ukb_genomicc_merged_fam[!, ["FID", "IID"]])
        end

        # Test report generation
        @test isfile(joinpath(results_dir, "call-make_report", "execution", "report.md"))
    end
end

end

true