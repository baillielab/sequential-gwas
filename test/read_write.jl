module TestReadWrite

using Test
using SequentialGWAS
using CSV
using DataFrames

PKGDIR = pkgdir(SequentialGWAS)
TESTDIR = joinpath(PKGDIR, "test")

@testset "Test complete-bim-with-ref" begin
    tmpdir = mktempdir()
    bim_file = joinpath(TESTDIR, "assets", "complete_bim_with_ref", "incomplete.bim")
    ref_bim_file = joinpath(TESTDIR, "assets", "qc_from_kgp", "kgp.bim")
    output = joinpath(tmpdir, "complete.bim")
    copy!(
        ARGS, [
            "complete-bim-with-ref", 
            bim_file, 
            ref_bim_file,
            "--output", output
    ]
    )
    julia_main()
    old_bim = SequentialGWAS.read_bim(bim_file)
    completed_bim = SequentialGWAS.read_bim(output)
    ref_bim = SequentialGWAS.read_bim(ref_bim_file)
    # Unchanged columns
    @test old_bim.CHR_CODE == completed_bim.CHR_CODE
    @test old_bim.BP_COORD == completed_bim.BP_COORD
    @test old_bim.POSITION == completed_bim.POSITION
    @test old_bim.ALLELE_2 == completed_bim.ALLELE_2
    # "." filled with missing allele
    @test all(completed_bim.ALLELE_1 .!== ".")
    # Update variant ID to chr:pos:ref:alt format
    @test all(x == 4 for x in length.(split.(completed_bim.VARIANT_ID, ":")))
    # Check a few ALLELE_1 updates
    variant_id = "chr1:56419369:C:T"
    ref_row = filter(:VARIANT_ID => ==(variant_id), ref_bim)
    completed_row = filter(:VARIANT_ID => ==(variant_id), completed_bim)
    @test ref_row.ALLELE_1 == completed_row.ALLELE_1
    variant_id = "chr2:22847824:G:A"
    ref_row = filter(:VARIANT_ID => ==(variant_id), ref_bim)
    completed_row = filter(:VARIANT_ID => ==(variant_id), completed_bim)
    @test ref_row.ALLELE_1 == completed_row.ALLELE_1
end

end

true