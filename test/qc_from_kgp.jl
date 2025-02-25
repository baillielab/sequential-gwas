module TestQCFilesFromKGP

using Test
using SequentialGWAS
using CSV
using DataFrames

TESTDIR = joinpath(pkgdir(SequentialGWAS), "test")

@testset "Test get_action" begin
    # After merging with KGP, if the variant is not in KGP, 
    # it will have a missing VARIANT_ID
    row = (KGP_VARIANT_ID = missing, )
    @test SequentialGWAS.get_action(row) == "DROP (VARIANT-NOT-IN-KGP)"
    # We only keep biallelic variants for now
    row = (KGP_VARIANT_ID = "chr1:1:A:C", NON_BIALLELIC = true, ) 
    @test SequentialGWAS.get_action(row) == "DROP (VARIANT-NON-BIALLELIC)"
    # If the reference allele is ALLELE_2 (major according to plink) and 
    # the alternative allele is ALLELE_2 (minor according to plink) then all is fine
    row = (KGP_VARIANT_ID = "chr1:1:A:C", NON_BIALLELIC = false, ALLELE_1 = "C", ALLELE_2 = "A")
    @test SequentialGWAS.get_action(row) == "KEEP (MINOR-MAJOR-MATCHING-ALT-REF)"
    row = (KGP_VARIANT_ID = "chr1:1:T:A", NON_BIALLELIC = false, ALLELE_1 = "A", ALLELE_2 = "T")
    @test SequentialGWAS.get_action(row) == "KEEP (MINOR-MAJOR-MATCHING-ALT-REF)"

    # There are two cases:
    ## (i) The variant is palindromic, then either:
    ##  - It wasn't properly flipped by GenomeStudio or earlier steps and needs to be flipped
    ##  - It has opposite allele frequency in our dataset (this should only happen when MAF ≈ 0.5 and we drop these cases)
    ## (ii) The variant is not palyndromic, it has an opposite allele frequency in our dataset (likely when MAF ≈ 0.5), we annotate it
    ## Non palindromic variant
    row = (KGP_VARIANT_ID = "chr1:1:C:A", NON_BIALLELIC = false, ALLELE_1 = "C", ALLELE_2 = "A")
    @test SequentialGWAS.get_action(row) == "KEEP (MINOR-MAJOR-REVERSED-ALT-REF)"
    ## Palindromic variants
    ### Balanced frequencies
    row = (KGP_VARIANT_ID = "chr1:1:T:A", NON_BIALLELIC = false, ALLELE_1 = "T", ALLELE_2 = "A", ALLELE_1_FREQ = 0.49, KGP_ALLELE_1_FREQ = 0.51)
    @test SequentialGWAS.get_action(row) == "DROP (PALINDROMIC-MINOR-MAJOR-REVERSED-ALT-REF-BALANCED)"
    ### Unbalanced frequencies, likely to be flipped
    row = (KGP_VARIANT_ID = "chr1:1:C:G", NON_BIALLELIC = false, ALLELE_1 = "C", ALLELE_2 = "G", ALLELE_1_FREQ = 0.3, KGP_ALLELE_1_FREQ = 0.7)
    @test SequentialGWAS.get_action(row) == "FLIP (PALINDROMIC-MINOR-MAJOR-REVERSED-ALT-REF)"

    # Otherwise, the alleles are non-matching the reference at all, 
    # this can be either because:
    # - The variant hasn't been flipped properly and we check the complementary alleles against the KGP 
    # (note again, the major/minor alleles may not correspond even once flipped and we report it)
    # - There is a problem with the reported alleles in the dataset and we drop the variant (this will also be reported)
    # These variants will be dropped
    row = (KGP_VARIANT_ID = "chr1:1:A:C", NON_BIALLELIC = false, ALLELE_1 = "T", ALLELE_2 = "A")
    @test SequentialGWAS.get_action(row) == "DROP (ALLELES-NOT-MATCHING-KGP)"
        row = (KGP_VARIANT_ID = "chr1:1:G:A", NON_BIALLELIC = false, ALLELE_1 = "G", ALLELE_2 = "C")
    @test SequentialGWAS.get_action(row) == "DROP (ALLELES-NOT-MATCHING-KGP)"
    # These variants will be flipped
    row = (KGP_VARIANT_ID = "chr1:1:A:C", NON_BIALLELIC = false, ALLELE_1 = "T", ALLELE_2 = "G")
    @test SequentialGWAS.get_action(row) == "FLIP (COMPLEMENT-MINOR-MAJOR-REVERSED-ALT-REF)"
    row = (KGP_VARIANT_ID = "chr1:1:T:G", NON_BIALLELIC = false, ALLELE_1 = "C", ALLELE_2 = "A")
    @test SequentialGWAS.get_action(row) == "FLIP (COMPLEMENT)"
end

end

true