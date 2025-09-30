module TestFineMapping

using GenomiccWorkflows
using Test
using CSV
using DataFrames

PKGDIR = pkgdir(GenomiccWorkflows)
TESTDIR = joinpath(PKGDIR, "test")

EXPECTED_FP_COLS = [
    "CHROM", "POS", "ID", "REF", "ALT", "LOCUS_ID", "PIP", "CS", "PHASED_R2"
]

@testset "Test make_clean_sample_file" begin
    tmpdir = mktempdir()
    # The file contains actual sample ids
    sample_file = joinpath(tmpdir, "sample_file.txt")
    open(sample_file, "w") do io
        write(io, "odap100\todap100")
    end
    @test GenomiccWorkflows.make_clean_sample_file(sample_file) == sample_file
    # The file contains a list of files containing sample ids written by REGENIE
    sample_file_1 = joinpath(tmpdir, "EUR.PHENO_1.chr1_PHENO_1.txt")
    write(sample_file_1, "odap1\todap1")
    sample_file_2 = joinpath(tmpdir, "EUR.PHENO_1.chr2_PHENO_1.txt")
    write(sample_file_2, "odap2\todap2")
    sample_file_3 = joinpath(tmpdir, "AFR.PHENO_1.chr1_PHENO_1.txt")
    write(sample_file_3, "odap3\todap3")
    sample_file_4 = joinpath(tmpdir, "AMR.PHENO_1.chr1_PHENO_1.txt")
    write(sample_file_4, "odap4\todap4")
    sample_file_5 = joinpath(tmpdir, "EUR.PHENO_2.chr1_PHENO_2.txt")
    write(sample_file_5, "odap5\todap5")

    sample_file = joinpath(tmpdir, "samples_file.txt")
    open(sample_file, "w") do io
        for filename in [sample_file_1, sample_file_2, sample_file_3, sample_file_4, sample_file_5]
            println(io, filename)
        end
    end

    output_sample_file = GenomiccWorkflows.make_clean_sample_file(sample_file; exclude=["AMR"], phenotype="PHENO_1", chr="1")
    @test readlines(output_sample_file) == ["odap1\todap1", "odap3\todap3"]
end


@testset "Test get_window_idxs" begin
    variant_ids = ["var0", "var1", "var2", "var3", "var4", "var5"]
    # LD variants have lower pos than lead: the lead is still included
    ld_variants = DataFrame(
        ID_A=["var3", "var3"],
        ID_B=["var1", "var2"]
    ) 
    idx_inf_sup = GenomiccWorkflows.get_window_idxs(
        ld_variants,
        variant_ids
    )
    @test idx_inf_sup == (2, 4)
    # LD variants have higher pos than lead: the lead is still included
    ld_variants = DataFrame(
        ID_A=["var3", "var3"],
        ID_B=["var4", "var5"]
    ) 
    idx_inf_sup = GenomiccWorkflows.get_window_idxs(
        ld_variants,
        variant_ids
    )
    @test idx_inf_sup == (4, 6)
    # LD variants are on both sides of the lead
    ld_variants = DataFrame(
        ID_A=["var3", "var3"],
        ID_B=["var2", "var4"]
    ) 
    idx_inf_sup = GenomiccWorkflows.get_window_idxs(
        ld_variants,
        variant_ids
    )
    @test idx_inf_sup == (3, 5)
end

@testset "Integration Test: finemap_significant_regions" begin
    gwas_results_file = joinpath(TESTDIR, "assets", "gwas", "results", "regenie.results.group.phenotype.tsv")
    pgen_prefix = joinpath(TESTDIR, "assets", "gwas", "imputed", "chr1.qced")
    sample_file = joinpath(TESTDIR, "assets", "gwas", "results", "sample_file.txt")
    covariates_file = joinpath(TESTDIR, "assets", "gwas", "covariates", "ukb_genomicc.covariates.csv")
    tmpdir = mktempdir()

    gwas_results = CSV.read(gwas_results_file, DataFrame)
    pvar = CSV.read(pgen_prefix * ".pvar", DataFrame; delim='\t', comment="##")
    # Identifying the variants that were not tested in GWAS
    GenomiccWorkflows.tag_variant_id_missing_from_gwas!(pvar, gwas_results)
    gwas_variants = filter(!endswith(".not_in_gwas"), pvar.ID)
    @test gwas_variants == [
        "chr1:14012312:T:C",
        "chr1:18100537:G:A",
        "chr1:22542609:T:C",
        "chr1:40310265:G:A",
        "chr1:92682820:C:T",
        "chr1:111622622:C:A"
    ]
    # Write those variants to a temporary PGEN fileset
    gwas_matched_pgen_prefix = joinpath(tmpdir, "chr1.gwas_matched")
    GenomiccWorkflows.write_new_pgen_from_gwas_results(pgen_prefix, gwas_matched_pgen_prefix, pvar, sample_file)
    @test countlines(gwas_matched_pgen_prefix * ".psam") == 1130
    @test countlines(gwas_matched_pgen_prefix * ".pvar") == 7
    # Write significant clumps: here none will be found but the function should still return an empty DataFrame
    output_prefix = joinpath(tmpdir, "fine_mapping_test")
    sig_clumps = GenomiccWorkflows.write_significant_clumps(gwas_matched_pgen_prefix, gwas_results_file;
        min_sig_clump_size = 1,
        output = string(output_prefix, ".clumps.tsv"),
        lead_pvalue = 1e-8,
        p2_pvalue = 1e-5,
        r2_threshold = 0.3,
        clump_kb = 100,
        clump_id_field = "ID",
        clump_pval_field = "LOG10P",
        allele_1_field = "ALLELE1"
    )
    @test nrow(sig_clumps) == 0
    @test names(sig_clumps) == [
        "#CHROM", "POS", "ID", "NEG_LOG10_P", "TOTAL", "NONSIG", 
        "S0.05", "S0.01", "S0.001", "S0.0001", "SP2"
    ]
    @test CSV.read(string(output_prefix, ".clumps.tsv"), DataFrame; delim="\t") == sig_clumps
    # Write significant clumps: here we set very permissive parameters to get a clump
    sig_clumps = GenomiccWorkflows.write_significant_clumps(gwas_matched_pgen_prefix, gwas_results_file;
        min_sig_clump_size = 1,
        output = string(output_prefix, ".clumps.tsv"),
        lead_pvalue = 0.3,
        p2_pvalue = 0.3,
        r2_threshold = 0.,
        clump_kb = 300_000,
        clump_id_field = "ID",
        clump_pval_field = "LOG10P",
        allele_1_field = "ALLELE1"
    )
    @test nrow(sig_clumps) == 1
    @test length(split(first(sig_clumps.SP2), ",")) == 2
    @test CSV.read(string(output_prefix, ".clumps.tsv"), DataFrame; delim="\t") == sig_clumps
    # To get two clumps
    sig_clumps = GenomiccWorkflows.write_significant_clumps(gwas_matched_pgen_prefix, gwas_results_file;
        min_sig_clump_size = 1,
        output = string(output_prefix, ".clumps.tsv"),
        lead_pvalue = 0.3,
        p2_pvalue = 0.3,
        r2_threshold = 0.,
        clump_kb = 70_000,
        clump_id_field = "ID",
        clump_pval_field = "LOG10P",
        allele_1_field = "ALLELE1"
    )
    # Test get_loci_to_finemap: no merging of clumps
    @test GenomiccWorkflows.get_loci_to_finemap(sig_clumps; window_kb=50) == [
        ["chr1:40310265:G:A", 0.88379, 40310265-50_000, 40310265+50_000],
        ["chr1:111622622:C:A", 0.935709, 111622622-50_000, 111622622+50_000]
    ]
    # Test get_loci_to_finemap: merging of clumps
    dist_between_clumps_kb = (111622622-40310265) / 1000
    @test GenomiccWorkflows.get_loci_to_finemap(sig_clumps; window_kb=dist_between_clumps_kb) == [
        ["chr1:111622622:C:A", 0.935709, 0.0, 1.82934979e8]
    ]
end

@testset "Integration Test: finemap_locus" begin
    finemap_window_kb = 30_000
    locus = ["chr1:40310265:G:A", 0.88379, 40310265-finemap_window_kb*1000, 40310265+finemap_window_kb*1000]
    locus_id, _, locus_start, locus_end = locus
    pgen_prefix = joinpath(TESTDIR, "assets", "gwas", "imputed", "chr1.qced")
    p = 3
    # Test function components
    ld_variants = GenomiccWorkflows.get_locus_variants_r2(locus_id, pgen_prefix; ld_window_kb=finemap_window_kb)
    @test nrow(ld_variants) == p
    X_df, variants_info = GenomiccWorkflows.dosages_from_pgen(pgen_prefix, ld_variants)
    @test Set(names(X_df)) == Set(vcat(ld_variants.ID_B, "chr1:40310265:G:A", "IID"))
    @test variants_info.ID == vcat(ld_variants.ID_B, "chr1:40310265:G:A")
    @test names(variants_info) == ["CHROM", "POS", "ID", "REF", "ALT"]
    for variant_id in variants_info.ID
        @test any(isnan.(X_df[!, 2])) == false
    end
    y_df = X_df[1:1000, [:IID]]
    y_df[!, :Y] = rand(Bool, nrow(y_df))
    X, y = GenomiccWorkflows.get_susie_inputs(X_df, y_df)
    @test size(X) == (1000, p+1)
    @test X isa Matrix{Float64}
    @test size(y) == (1000,)
    @test y isa Vector{Float64}
    finemapping_results = GenomiccWorkflows.susie_finemap(X, y; n_causal=2)
    @test length(finemapping_results[:pip]) == 4
    # Post processing
    GenomiccWorkflows.postprocess_finemapping_results!(variants_info, finemapping_results, ld_variants)
    @test names(variants_info) == EXPECTED_FP_COLS
    @test nrow(variants_info) == p+1
    @test any(ismissing.(variants_info.PHASED_R2)) == false
    # Test full function
    finemapping_results = GenomiccWorkflows.finemap_locus(locus, pgen_prefix, y_df;
        Xtype="dosages",
        n_causal=2,
        finemap_window_kb=finemap_window_kb,
    )
    @test names(finemapping_results) == EXPECTED_FP_COLS
    @test nrow(finemapping_results) == p+1
    # Test genotypes_from_pgen
    X_df, variants_info = GenomiccWorkflows.genotypes_from_pgen(pgen_prefix, ld_variants)
    @test Set(names(X_df)) == Set(vcat(ld_variants.ID_B, "chr1:40310265:G:A", "IID"))
    @test variants_info.ID == vcat(ld_variants.ID_B, "chr1:40310265:G:A")
    @test names(variants_info) == ["CHROM", "POS", "ID", "REF", "ALT"]
    for variant_id in variants_info.ID
        @test any(isnan.(X_df[!, 2])) == false
    end
end

@testset "Integration Test: finemap_locus_rss" begin
    finemap_window_kb = 30_000
    locus = ["chr1:40310265:G:A", 0.88379, 40310265-finemap_window_kb*1000, 40310265+finemap_window_kb*1000]
    locus_id, _, locus_start, locus_end = locus
    pgen_prefix = joinpath(TESTDIR, "assets", "gwas", "imputed", "chr1.qced")
    p = 3
    n = 10
    # Get LD matrix
    locus_id, _ = locus
    ld_variants = GenomiccWorkflows.get_locus_variants_r2(locus_id, pgen_prefix; ld_window_kb=finemap_window_kb)
    R, variants = GenomiccWorkflows.get_LD_matrix(pgen_prefix, locus)
    @test R isa Matrix{Float64}
    @test size(R) == (p+1, p+1)
    @test variants == ["chr1:14012312:T:C", "chr1:18100537:G:A", "chr1:22542609:T:C", "chr1:40310265:G:A"]
    # Run susie_rss_finemap
    gwas_results = DataFrame(
        CHROM = ["1", "1", "1", "1", "1"],
        ID = ["chr1:14012312:T:C", "chr1:18100537:G:A", "chr1:22542609:T:C", "chr1:40310265:G:A", "toto"],
        GENPOS = [14012312, 18100537, 22542609, 40310265, 111622622],
        ALLELE0 = ["T", "G", "T", "G", "C"],
        ALLELE1 = ["C", "A", "C", "A", "A"],
        BETA = randn(5),
        SE = rand(5),
    )
    variants_info = GenomiccWorkflows.initialize_variants_info_rss(pgen_prefix, variants, gwas_results)
    @test variants_info.ID == ["chr1:14012312:T:C", "chr1:18100537:G:A", "chr1:22542609:T:C", "chr1:40310265:G:A"]
    @test names(variants_info) == ["CHROM", "POS", "ID", "REF", "ALT", "BETA", "SE"]
    # Fake phenotype vector
    y = rand(n)
    finemapping_results = GenomiccWorkflows.susie_rss_finemap(R, variants_info, y; n_causal=1)
    @test length(finemapping_results[:pip]) == 4
    # Post processing
    GenomiccWorkflows.postprocess_finemapping_results!(variants_info, finemapping_results, ld_variants)
    @test names(variants_info) == EXPECTED_FP_COLS
    @test any(ismissing.(variants_info.PHASED_R2)) == false
    # Full function
    finemapping_results = GenomiccWorkflows.finemap_locus_rss(locus, gwas_results, pgen_prefix, y;
        n_causal=1,
        finemap_window_kb=finemap_window_kb,
    )
    @test names(finemapping_results) == EXPECTED_FP_COLS
    @test nrow(finemapping_results) == p+1
end

@testset "Test region_plot" begin
    finemapping_results = GenomiccWorkflows.harmonize_finemapping_results(
            CSV.read(
        joinpath(TESTDIR, "assets", "gwas", "results", "results.all_chr.EUR.SEVERE_COVID_19.finemapping.tsv"), 
        DataFrame; 
        delim="\t"
    ))
    finemapping_results = finemapping_results[finemapping_results.LOCUS_ID .== "rs7515509", :]
    gwas_results = GenomiccWorkflows.harmonize_gwas_results(
            CSV.read(
        joinpath(TESTDIR, "assets", "gwas", "results", "results.all_chr.EUR.SEVERE_COVID_19.gwas.tsv"), 
        DataFrame; 
        delim="\t"
    ))
    region_data = innerjoin(
        gwas_results,
        DataFrames.select(finemapping_results, [:ID, :REF, :ALT, :PIP, :CS, :LOCUS_ID, :PHASED_R2]), 
        on=[:ID]
    )
    fig = GenomiccWorkflows.region_plot(region_data)
    @test fig !== nothing
end


end

true