using SequentialGWAS #hide
using CSV #hide
using DataFrames #hide
using MarkdownTables #hide
using Markdown #hide
using CairoMakie #hide
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
function qc_report(log_file) #hide
    filtered = Dict() #hide
    for line in readlines(log_file) #hide
        if startswith(line, "--geno") && occursin("removed due to missing genotype data", line) #hide
            filtered["--geno"] = parse(Int, split(line, " ")[2]) #hide
        elseif startswith(line, "--hwe") && occursin("removed due to Hardy-Weinberg exact test", line) #hide
            filtered["--hwe"] = parse(Int, split(line, " ")[2]) #hide
        elseif startswith(line, "--rm-dup") && occursin("duplicated", line) #hide
            filtered["--rm-dup"] = parse(Int, split(line, " ")[2]) #hide
        elseif occursin("removed due to missing genotype data (--mind)", line) #hide
            filtered["--mind"] = parse(Int, split(line, " ")[1]) #hide
        end #hide
    end #hide
    return [get(filtered, "--geno", 0), get(filtered, "--hwe", 0), get(filtered, "--rm-dup", 0), get(filtered, "--mind", 0)] #hide
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
    "" => ["--geno", "--hwe", "--rm-dup", "--mind"], #hide
    "r8 release" => qc_report(release_r8_qc_logs), #hide
    "2021 - 2023 release" => qc_report(release_2021_2023_qc_logs), #hide
    "2024 - now release " => qc_report(release_2024_qc_logs) #hide
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

Samples dropped from each release due to overlap with a higher priority release.
=#

println("r8 release: ", countlines(dup_samples_r8)) #hide
println("2021 - 2023 release: ", countlines(dup_samples_2021_2023)) #hide
println("2024 - now release: ", countlines(dup_samples_2024_now)) #hide

#=
## Whole Genome Sequencing (Optional)

The whole-genome sequencing data is then genotyped to include the variants present in the genotyping arrays.

Number of variants and individuals in filesets.
=#

function wgs_table(wgs_prefix) #hide
    if wgs_prefix !== "NO_WGS_SAMPLES" #hide
        table = DataFrame(Dict( #hide
            "" => ["Variants", "Individuals"], #hide
            "WGS" => count_variants_and_indv(wgs_prefix), #hide
        )) #hide
        return table |> markdown_table() #hide
    else #hide
        return nothing #hide
    end #hide
end #hide

wgs_table(wgs_prefix) #hide

#=
## Merged Datasets

We then proceed to merging the genotyping arrays and optionally WGS datasets.

### After Merge

Number of variants and individuals in fileset after merge.
=#

table = DataFrame(Dict( #hide
    "" => ["Variants", "Individuals"], #hide
    "Merged" => count_variants_and_indv(merged_genotypes_prefix), #hide
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
    "MergedQCed" => count_variants_and_indv(merged_qced_genotypes_prefix), #hide
)) #hide
table |> markdown_table() #hide

#=
### After PCA QC (Final Step)

Some variants may have high PCA loadings, we remove them.

- Loadings before Removal
=#

cp(string(pca_plot_prefix, ".before_pca_qc.loadings.png"), "loadings_before_pca.png") #hide
MD(string("[loadings-before](loadings_before_pca.png)")) #hide

# - PCA / Ancestry estimates before Removal
cp(string(pca_plot_prefix, ".before_pca_qc.all.png"), "pca_ancestry_before.png") #hide
MD(string("[pca-ancestry-before](pca_ancestry_before.png)")) #hide

# - Number of variants excluded by PCA QC (high loadings)

countlines(high_loadings_variants) #hide

# - Loadings after removal
cp(string(pca_plot_prefix, ".after_pca_qc.loadings.png"), "loadings_after_pca.png") #hide
MD(string("[loadings-after](loadings_after_pca.png)")) #hide

# - PCA / Ancestry estimates
cp(string(pca_plot_prefix, ".after_pca_qc.all.png"), "pca_ancestry_after.png") #hide
MD(string("[pca-ancestry-after](pca_ancestry_after.png)")) #hide

#=
## Some Statistics

- Number of variants and individuals in final fileset

These are the final number of variants and individuals in the merged dataset.
=#

table = DataFrame(Dict( #hide
    "" => ["Variants", "Individuals"], #hide
    "Final" => count_variants_and_indv(final_genotypes_prefix), #hide
)) #hide
table |> markdown_table() #hide

#=
- Estimated Ancestries

Number of individuals per estimated ancestry group
=#

covariates = CSV.read(covariates_file, DataFrame) #hide
ancestries = combine(groupby(covariates, :ANCESTRY), nrow => :N) #hide
ancestries |> markdown_table() #hide

#=
- Samples per platform

Number of individuals per genetic measurement platform
=#

platforms = combine(groupby(covariates, :PLATFORM), nrow => :N) #hide
platforms |> markdown_table() #hide