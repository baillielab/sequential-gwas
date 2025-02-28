
using SequentialGWAS #hide
using CSV #hide
using DataFrames #hide
using MarkdownTables #hide
using Markdown #hide
#hide
struct MD str end                                                      #hide
Base.show(io::IO, ::MIME"text/markdown", md::MD) = print(io, md.str)   #hide
function count_variants_and_indv(file_prefix) #hide
    bim_file = string(file_prefix, ".bim") #hide
    fam_file = string(file_prefix, ".fam") #hide
    return [countlines(bim_file), countlines(fam_file)] #hide
end #hide
function count_dropped_variants_and_indiv(file_prefix) #hide
    variants = CSV.read(string(file_prefix, ".filtered_variants.csv"), DataFrame) #hide
    individuals = CSV.read(string(file_prefix, ".filtered_samples.csv"), DataFrame) #hide
    return [nrow(variants), nrow(individuals)] #hide
end #hide
#hide
nothing #hide
#=
# Genomicc Aggregation Pipeline Report

## Array Genotypes

This section presents the effect of each step applied to genotyping arrays in the pipeline.

### LiftOver

Only the release r8 and the 2021-2023 release are in GRCh37 build. They are lifted over, which can result in a loss of
genetic variants (those that cannot be mapped to GRCh38). 

Variants dropped by liftover procedure.

- Release r8
=#

countlines(unlifted_r8) #hide

# - Release 2021 - 2023

countlines(unlifted_2021_2023) #hide

#=

### Bi-allelic only

We only keep bi-allelic variants at the moment. This is the total number of variants and individuals in each fileset.
=#

table = DataFrame(Dict( #hide
    "" => ["Variants", "Individuals"], #hide
    "r8 release" => count_variants_and_indv(initial_bed_prefix_r8), #hide
    "2021 - 2023 release" => count_variants_and_indv(initial_bed_prefix_2021_2023), #hide
    "2024 - now release" => count_variants_and_indv(initial_bed_prefix_2024_now) #hide
)) #hide
table |> markdown_table() #hide

#=

### Basic QC

We apply basic QC filters to drop variants and individuals with poor genotyping. 

Variants and individuals dropped by basic QC (plink2 --mind --geno --hwe).
=#

table = DataFrame(Dict( #hide
    "" => ["Variants", "Individuals"], #hide
    "r8 release" => count_dropped_variants_and_indiv(basic_qc_prefix_r8), #hide
    "2021 - 2023 release" => count_dropped_variants_and_indiv(basic_qc_prefix_2021_2023), #hide
    "2024 - now release " => count_dropped_variants_and_indiv(basic_qc_prefix_2024_now) #hide
)) #hide
table |> markdown_table() #hide

#=
### Shared variants and individuals across releases and WGS 

Only variants present in all releases are kept.

- Shared variants across releases.
=#

countlines(shared_variants) #hide

#=
- Shared individuals across releases.

Some individuals may be present in multiple releases (and have whole-genome sequencing). We give priority to data sources in the following order: WGS, 2024 - now, 2021 - 2023, r8 release.

Variants dropped from each release due to overlap with a higher priority release.
=#

println("r8 release: ", countlines(dup_samples_r8)) #hide
println("2021 - 2023 release: ", countlines(dup_samples_2021_2023)) #hide
println("2024 - now release: ", countlines(dup_samples_2024_now)) #hide

#=
## Whole Genome Sequencing

The whole-genome sequencing data is then genotyped to include the variants present in the genotyping arrays.

Number of variants and individuals in filesets.
=#

table = DataFrame(Dict( #hide
    "" => ["Variants", "Individuals"], #hide
    "WGS" => count_variants_and_indv(wgs_prefix), #hide
)) #hide
table |> markdown_table() #hide

#=
## Merged Datasets

We then proceed to merging the geotyping array and WGS datasets.

### After Merge

Number of variants and individuals in fileset after merge.
=#

table = DataFrame(Dict( #hide
    "" => ["Variants", "Individuals"], #hide
    "WGS" => count_variants_and_indv(merged_genotypes_prefix), #hide
)) #hide
table |> markdown_table() #hide

#=
### Relatedness

Some individuals may be strongly related biasing downstream analyses. We will remove them.

Number of unrelated individuals
=#

countlines(unrelated_individuals) #hide

#=
### After Basic QC

We re-apply basic QC filters (plink2 --mind --geno --hwe) and only keep unrelated individuals 

Number of variants and individuals in merged fileset after QC.
=#

table = DataFrame(Dict( #hide
    "" => ["Variants", "Individuals"], #hide
    "WGS" => count_variants_and_indv(merged_qced_genotypes_prefix), #hide
)) #hide
table |> markdown_table() #hide

#=
### After PCA QC

Some variants may have high PCA loadings, we remove them.

- Loadings before Removal
=#

cp(string(pca_plot_prefix, ".before_pca_qc.loadings.png"), "loadings_before_pca.png") #hide
MD(string("[loadings-before](loadings_before_pca.png)")) #hide

# - Number of variants excluded by PCA QC (high loadings)

countlines(high_loadings_variants) #hide

# - Loadings after removal
cp(string(pca_plot_prefix, ".after_pca_qc.loadings.png"), "loadings_after_pca.png") #hide
MD(string("[loadings-after](loadings_after_pca.png)")) #hide

# - PCA / Ancestry estimates
cp(string(pca_plot_prefix, ".after_pca_qc.all.png"), "pca_ancestry.png") #hide
MD(string("[pca-ancestry](pca_ancestry.png)")) #hide

#=
- Number of variants and individuals in final fileset

These are the final number of variants and individuals in the merged dataset.
=#

table = DataFrame(Dict( #hide
    "" => ["Variants", "Individuals"], #hide
    "WGS" => count_variants_and_indv(final_genotypes_prefix), #hide
)) #hide
table |> markdown_table() #hide


