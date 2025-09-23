module TestGWASPLots
using Test
using GenomiccWorkflows

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