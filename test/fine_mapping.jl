module TestFineMapping

using GenomiccWorkflows
using Test
using CSV
using DataFrames

PKGDIR = pkgdir(GenomiccWorkflows)
TESTDIR = joinpath(PKGDIR, "test")

@testset "Integration Test" begin
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
    @test CSV.read(string(output_prefix, ".clumps.tsv"), DataFrame; delim=",") == sig_clumps
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
    @test CSV.read(string(output_prefix, ".clumps.tsv"), DataFrame; delim=",") == sig_clumps
    # Get variants in LD with clump lead
    sample_list = getindex.(split.(readlines(sample_file), "\t"), 2)
    y = GenomiccWorkflows.get_phenotype(covariates_file, sample_list, "SEVERE_COVID_19")
    clump = first(sig_clumps)

    clump_finemapping_results = GenomiccWorkflows.finemap_clump(clump.ID, gwas_matched_pgen_prefix, y, sample_list;
        n_causal=1,
        ld_window_kb=300_000,
        ld_window_r2=0,
    )
    @test names(clump_finemapping_results) == [
        "#CHROM", "POS", "ID", "REF", "ALT", "PIP", "CLUMP_ID", "CS"
    ]
    @test nrow(clump_finemapping_results) == 5
end


end

true