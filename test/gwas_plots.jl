module TestGWASPLots
using Test
using GenomiccWorkflows
using DataFrames

@testset "Test harmonize_gwas_results" begin
    # If the LOG10P contains NA it will be read as a string column
    ## The NaNs are filtered for plotting
    results = DataFrame(
        LOG10P = ["NA", "1", "1", "2"],
        CHROM = [1, 1, 2, 2],
        GENPOS = [1000, 2000, 3000, 4000],
        ID = ["rs1", "rs2", "rs3", "rs4"],
        A1FREQ = ["0.1", "0.2", "NA", "0.4"]
    )
    harmonized_results = GenomiccWorkflows.harmonize_gwas_results(results)
    @test harmonized_results.CHR == ["1", "1", "2", "2"]
    @test harmonized_results.BP == results.GENPOS
    @test harmonized_results.SNP == results.ID
    @test harmonized_results.P[1] === NaN
    @test harmonized_results.P[2:end] == [0.1, 0.1, 0.01]

    # If the LOG10P has no NA it will be read as a float column
    results = DataFrame(
        LOG10P = [1, 1, 2, 2],
        CHROM = [1, 1, 2, 2],
        GENPOS = [1000, 2000, 3000, 4000],
        ID = ["rs1", "rs2", "rs3", "rs4"],
        A1FREQ = [0.1, 0.2, 0.3, 0.4]
    )
    harmonized_results = GenomiccWorkflows.harmonize_gwas_results(results)
    @test harmonized_results.CHR == ["1", "1", "2", "2"]
    @test harmonized_results.BP == results.GENPOS
    @test harmonized_results.SNP == results.ID
    @test harmonized_results.P == [0.1, 0.1, 0.01, 0.01]
end

@testset "Test gwas_plots" begin
    tmpdir = mktempdir()
    output_prefix = joinpath(tmpdir, "plot")
    copy!(ARGS,[
        "gwas-plots", 
        "test/assets/gwas/results/results.all_chr.EUR.SEVERE_COVID_19.gwas.tsv", 
        "test/assets/gwas/results/results.all_chr.EUR.SEVERE_COVID_19.finemapping.tsv", 
        "--maf=0.01",
        "--output-prefix=$output_prefix"
        ]
    )
    julia_main()
    @test isfile(string(output_prefix, ".EUR.SEVERE_COVID_19.manhattan.png"))
    @test isfile(string(output_prefix, ".EUR.SEVERE_COVID_19.qq.png"))
    @test isfile(string(output_prefix, ".EUR.SEVERE_COVID_19.rs12732514.locuszoom.png"))
    @test isfile(string(output_prefix, ".EUR.SEVERE_COVID_19.rs7515509.locuszoom.png"))
end

end

true