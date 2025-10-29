#=

This script can be used to generate a subset of the 1000GP. 

- Variants are filtered to match those in the `GRC38_map_file`.
- Individuals are filtered to keep `P` individuals within each ancestry.

Arguments:

- GP_DIR: Where the 1000GP files are located
- MAP_FILE: The map file to use for variants filtering
- OUTDIR: Where the files will be output
- P: Number of 1000GP individuals to keep per ancestry/gender
=#
using GenomiccWorkflows
using DataFrames
using CSV

GP_DIR = joinpath("assets", "resources", "thousandsGP")
OUT_DIR = joinpath("test", "assets", "kgp")
MAP_FILE_GENOMIC = joinpath("test", "assets", "genomicc", "genotyping_arrays", "mock.release_2024_now.map")
P = 20

# Write individuals subset
panel_file_1000_gp = CSV.read(
    joinpath(GP_DIR, "20130606_g1k_3202_samples_ped_population.txt"), 
    DataFrame, 
    ntasks=1
)

subpop_file = joinpath(OUT_DIR, "individuals_subset.csv")
subpop = mapreduce(vcat, groupby(panel_file_1000_gp, [:Superpopulation, :Sex])) do group
    DataFrames.select(group[1:P, :], :SampleID)
end
CSV.write(subpop_file, subpop, delim="\t", header=false)

# Write VCF subsets

map_file_genomic = DataFrame(read_map(MAP_FILE_GENOMIC), [:CHR, :ID, :POS, :BP])
map_file_genomic.CHR .= "chr" .* string.(map_file_genomic.CHR)
map_file_genomic.BP_START = map_file_genomic.BP .- 10
map_file_genomic.BP_END = map_file_genomic.BP .+ 10
(chr_key, chr_group) = first(pairs(groupby(map_file_genomic, :CHR)))
for (chr_key, chr_group) in pairs(groupby(map_file_genomic, :CHR))
    chr = chr_key.CHR
    @info(string("Processing chromosome ", chr))
    # Creating variables
    filename = string("CCDG_14151_B01_GRM_WGS_2020-08-05_", chr, ".filtered.shapeit2-duohmm-phased.vcf.gz")
    input_vcf = joinpath(GP_DIR, filename)
    @assert isfile(input_vcf)
    output_vcf = joinpath(OUT_DIR, filename)
    # Write variants to filter to file
    variants_file = joinpath(OUT_DIR, string("variants_", chr, ".txt"))
    CSV.write(
        variants_file, 
        DataFrames.select(chr_group, [:CHR, :BP_START, :BP_END]), 
        delim="\t", 
        header=false
    )
    # Filter and Index
    run(`bcftools view -Oz -R $variants_file $input_vcf -S $subpop_file -o $output_vcf`)
    run(`tabix -p vcf $output_vcf`)
    # Clean
    rm(variants_file)
end

rm(subpop_file)
