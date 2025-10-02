module TestMetaAnalysis

using Test
using GenomiccWorkflows
using DataFrames
using CSV

PKGDIR = pkgdir(GenomiccWorkflows)
TESTDIR = joinpath(PKGDIR, "test")

@testset "Test load_meta_analysis_worklist" begin
    tmpdir = mktempdir()
    gwas_results_list_file = joinpath(tmpdir, "regenie_files_list.txt")
    regenie_files_list = readdir(joinpath(TESTDIR, "assets", "gwas", "meta_analysis"), join=true)
    open(gwas_results_list_file, "w") do io
        for f in regenie_files_list
            println(io, f)
        end
    end
    # No exclusion of groups
    exclude = []
    worklist = GenomiccWorkflows.load_meta_analysis_worklist(gwas_results_list_file; exclude = exclude)
    @test sort(worklist[!, [:PHENOTYPE, :GROUP]]) == DataFrame(
        PHENOTYPE = ["SEVERE_COVID_19", "SEVERE_PNEUMONIA", "SEVERE_PNEUMONIA", "SEVERE_PNEUMONIA"],
        GROUP = ["EAS", "AFR", "EUR", "SAS"]
    )
    # No exclusion of groups
    exclude =["EAS", "AFR"]
    worklist = GenomiccWorkflows.load_meta_analysis_worklist(gwas_results_list_file; exclude = exclude)
    @test sort(worklist[!, [:PHENOTYPE, :GROUP]]) == DataFrame(
        PHENOTYPE = ["SEVERE_PNEUMONIA", "SEVERE_PNEUMONIA"],
        GROUP = ["EUR", "SAS"]
    )
end

@testset "Test append_GWAS_info_to_meta_analysis_results!" begin
    phenotype_gwas_files = filter(
        x -> occursin("SEVERE_PNEUMONIA", x),
        readdir(joinpath(TESTDIR, "assets", "gwas", "meta_analysis"), join=true)
    )
    metal_results = CSV.read(
        joinpath(TESTDIR, "assets", "gwas", "meta_analysis", "AFR.SEVERE_PNEUMONIA.gwas.tsv"), 
        DataFrame, 
        delim="\t", 
        select=[:ID]
    )
    GenomiccWorkflows.append_GWAS_info_to_meta_analysis_results!(metal_results, phenotype_gwas_files)
    @test nrow(metal_results) == 19
    @test names(metal_results) == [
        "ID", "CHROM", "GENPOS", "ALLELE0", "ALLELE1", "A1FREQ", "N", "NGROUPS"
    ]
    # "chr1:14012312:T:C" is not in the EUR file
    @test metal_results[metal_results.ID .== "chr1:14012312:T:C", :NGROUPS] == [2]
    
    for gwas_file in phenotype_gwas_files
        gwas_results = CSV.read(gwas_file, DataFrame; delim="\t", select=[:ID, :CHROM, :GENPOS, :ALLELE0, :ALLELE1, :A1FREQ, :N])
        joined = leftjoin(
            select(metal_results, :ID, :A1FREQ => :MINA1FREQ, :N => :SUM_N), 
            select(gwas_results, :ID, :A1FREQ => :GROUP_A1FREQ, :N => :GROUP_N), 
            on=:ID
        )
        @test all(skipmissing(joined.MINA1FREQ .<= joined.GROUP_A1FREQ))
        @test all(skipmissing(joined.SUM_N .>= joined.GROUP_N))
    end
end

@testset "Test meta_analyse" begin
    tmpdir = mktempdir()
    output_prefix = joinpath(tmpdir, "gwas.meta_analysis")
    gwas_results_list_file = joinpath(tmpdir, "regenie_files_list.txt")
    regenie_files_list = readdir(joinpath(TESTDIR, "assets", "gwas", "meta_analysis"), join=true)
    open(gwas_results_list_file, "w") do io
        for f in regenie_files_list
            println(io, f)
        end
    end
    copy!(ARGS, 
        ["meta-analyse", 
        gwas_results_list_file,
        "--exclude=SAS",
        "--output-prefix=$output_prefix"
    ])
    julia_main()
    expected_cols = Set([
        "ID", "BETA", "SE", "LOG10P", "DIRECTION", 
        "HET_ISQ", "HET_CHISQ", "HET_DF", "LOG10P_HET", 
        "CHROM", "GENPOS", "ALLELE0", "ALLELE1", 
        "A1FREQ", "N", "NGROUPS"
    ])
    meta_covid_19 = CSV.read(joinpath(tmpdir, "gwas.meta_analysis.SEVERE_COVID_19.gwas.tsv"), DataFrame)
    @test Set(names(meta_covid_19)) == expected_cols
    @test all(meta_covid_19.NGROUPS .== 1)
    
    meta_pneumonia = CSV.read(joinpath(tmpdir, "gwas.meta_analysis.SEVERE_PNEUMONIA.gwas.tsv"), DataFrame)
    @test Set(names(meta_pneumonia)) == expected_cols
    meta_pneumonia[meta_pneumonia.ID .== "chr1:14012312:T:C", :NGROUPS] == [1]
    @test all(meta_pneumonia.NGROUPS .<= 2)
end

end

true