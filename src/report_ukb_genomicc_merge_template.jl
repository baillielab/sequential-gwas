using GenomiccWorkflows #hide
using DataFrames #hide
using MarkdownTables #hide
using Markdown #hide
using CSV #hide
function find_header_line(file) #hide
    open(file) do io #hide
        line_idx = 1 #hide
        while true #hide
            line = readline(io) #hide
            if startswith(line, "#CHROM") #hide
                return line_idx #hide
            end #hide
            line_idx += 1 #hide
        end #hide
    end #hide
end #hide
function make_imputed_counts_table(ukb_genomicc_imputed_files_list) #hide
    imputed_df = DataFrame(FILE=readlines(ukb_genomicc_imputed_files_list)) #hide
    imputed_df.CHR = map(x -> split(x, ".")[end-1], imputed_df.FILE) #hide
    chr_counts = [] #hide
    for group in groupby(imputed_df, :CHR) #hide
        @assert nrow(group) == 2 #hide
        psam = only(filter(x -> endswith(x, ".psam"), group.FILE)) #hide
        n_samples = countlines(psam) - 1 #hide
        pvar = only(filter(x -> endswith(x, ".pvar"), group.FILE)) #hide
        n_variants = countlines(pvar) - find_header_line(pvar) #hide
        push!(chr_counts, (CHR=group.CHR[1], N_VARIANTS=n_variants, N_SAMPLES=n_samples)) #hide
    end #hide
    return DataFrame(chr_counts) |> markdown_table() #hide
end #hide
function ancestry_counts(ukb_genomicc_covariates_file) #hide
    covariates_df = CSV.read(ukb_genomicc_covariates_file, DataFrame) #hide
    counts = combine( #hide
        groupby(covariates_df, [:ANCESTRY_ESTIMATE]), #hide
        :COHORT => (col -> sum(x == "GENOMICC" for x in col)) => :GENOMICC, #hide
        :COHORT => (col -> sum(x == "UKB" for x in col)) => :UKB, #hide
        nrow => :TOTAL #hide
    ) #hide
    return counts |> markdown_table() #hide
end #hide
nothing #hide
#=
# UK Biobank & Genomicc Merge Workflow Report

The GenOMICC cohort is concatenated with the UK Biobank cohort to create a case/control dataset suitable for genetic analyses of critical care susceptibility. 
Only UK Biobank individuals with no record of being admitted to critical care [see this](https://biobank.ndph.ox.ac.uk/ukb/rectab.cgi?id=1065) are retained in the merged dataset.
This report summarizes the merging process, including the number of variants and individuals retained after critical steps.

## Genotypes

Because the UK Biobank and GenOMICC genotyping arrays differ, we instead filter the UK Biobnak imputed genotypes to match variants in the GenOMICC dataset.

Number of variants after merge:
=#
countlines(ukb_genomicc_merged_bim_file) #hide
#=
Number of individuals after merge:
=#
countlines(ukb_genomicc_merged_fam_file) #hide
#=
## Imputed Genotypes

The UK Biobank and GenOMICC TOPMed imputed genotypes are merged together, still filtering UK Biobank individuals admitted to critical care.
=#
make_imputed_counts_table(ukb_genomicc_imputed_files_list) #hide

#=
## Covariates

The UK Biobank and GenOMICC covariates are merged together, still filtering UK Biobank individuals admitted to critical care.

Repartition of individuals per ancestry
=#
ancestry_counts(ukb_genomicc_covariates_file) #hide

#=
Total number of individuals
=#
countlines(ukb_genomicc_covariates_file) - 1 #hide