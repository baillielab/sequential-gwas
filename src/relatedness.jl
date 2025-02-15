"""

Only keeps the first individual within each family and writes them to `outfile`.
"""
function kgp_unrelated_individuals(pedigree_file; outfile="kgp_samples_to_keep.txt")
    pedigrees = CSV.read(pedigree_file, DataFrame)
    unrelated_individuals = mapreduce(
        group -> group[[1], [:FamilyID, :SampleID]], 
        vcat, 
        groupby(pedigrees, :FamilyID)
    )
    CSV.write(outfile, unrelated_individuals[!, [:SampleID]], delim='\t', header=false)
end