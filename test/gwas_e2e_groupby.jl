module TestGWASE2E2

# This end to end test runs with the following conditions:
# - a grouped analysis defined by the SUPERPOPULATION column
# - two binary phenotypes: SEVERE_COVID_19 and SEVERE_PNEUMONIA
# - finemapping threhsolds have been lowered to be effectively performed
# - meta analysis is performed and the AMR group is excluded from it

using Test
using GenomiccWorkflows
using DataFrames
using CSV

PKGDIR = pkgdir(GenomiccWorkflows)
TESTDIR = joinpath(PKGDIR, "test")

function dir_contains_subdir(dir_name, subdir_name)
    subdirs = readdir(dir_name)
    if subdir_name in subdirs
        return true
    else
        first_dir = joinpath(dir_name, first(subdirs))
        if isdir(first_dir)
            return dir_contains_subdir(first_dir, subdir_name)
        else 
            return false
        end
    end
end

cmd_args = haskey(ENV, "CROMWELL_PATH") ?
    ["-jar", ENV["CROMWELL_PATH"]] :
    ["-Dconfig.file=conf/cromwell.mac.conf", "-jar", "/Users/olabayle/cromwell/cromwell-90.jar"]

cmd = Cmd([
    "java", cmd_args...,
    "run", "rap_workflows/gwas/workflow.wdl",
    "--inputs", joinpath(TESTDIR, "assets", "gwas.bygroup.json"),
])
rc = run(cmd)
@test rc.exitcode == 0

results_dirs = readdir("cromwell-executions/gwas", join=true)
results_dir = results_dirs[argmax(mtime(d) for d in results_dirs)]
expected_groups = Set([
    "AFR.SEVERE_PNEUMONIA", 
    "AMR.SEVERE_PNEUMONIA", 
    "EAS.SEVERE_COVID_19", 
    "EAS.SEVERE_PNEUMONIA", 
    "EUR.SEVERE_PNEUMONIA",
    "SAS.SEVERE_PNEUMONIA"
])
# Test groups and covariates: first see which groups make the case/control constraint
groups_prep_dir = joinpath(results_dir, "call-make_covariates_and_groups", "execution")
covariates = CSV.read(joinpath(groups_prep_dir, "gwas.covariates.csv"), DataFrame)
@test "AGE_x_AGE" in names(covariates)
covid_19_groups_not_passing_cc_threshold = Set(filter(
    x -> x.nrow < 600,
    sort(combine(groupby(covariates, [:SUPERPOPULATION, :SEVERE_COVID_19], skipmissing=true), nrow), :nrow)
).SUPERPOPULATION)
@test covid_19_groups_not_passing_cc_threshold == Set(["ADMIXED", "EUR", "AMR", "AFR", "SAS"])
pneumonia_groups_not_passing_cc_threshold = Set(filter(
    x -> x.nrow < 600,
    sort(combine(groupby(covariates, [:SUPERPOPULATION, :SEVERE_PNEUMONIA], skipmissing=true), nrow), :nrow)
).SUPERPOPULATION)
@test pneumonia_groups_not_passing_cc_threshold == Set(["ADMIXED"])
# Now check the groups files
for (ancestry, phenotype) in Iterators.product(
        ["AFR", "AMR", "EAS", "EUR", "SAS"],
        ["SEVERE_COVID_19", "SEVERE_PNEUMONIA"]
    )
    potential_sample_list = joinpath(groups_prep_dir, "gwas.individuals.$ancestry.$phenotype.txt")
    if "$ancestry.$phenotype" in expected_groups
        @test isfile(potential_sample_list)
    else
        @test !isfile(potential_sample_list)
    end
end
covariate_list = readlines(joinpath(groups_prep_dir, "gwas.covariates_list.txt"))
@test Set(covariate_list) == Set(["AGE", "SEX", "AGE_x_AGE"])

# Test BED groups qced
bed_dir = joinpath(results_dir, "call-make_group_bed_qced")
groups = Set([])
for shard in [0, 1, 2, 3, 4, 5] # expected 6 shards for 6 sample lists
    execution_dir = joinpath(bed_dir, "shard-$shard", "execution")
    files = readdir(execution_dir)
    fam_file = files[findfirst(endswith(".fam"), files)]
    push!(groups, splitext(fam_file)[1])
end
@test groups == expected_groups

# Test LD pruning
ld_prune_dir = joinpath(results_dir, "call-groups_ld_prune")
groups = Set([])
for shard in [0, 1, 2, 3, 4, 5]
    execution_dir = joinpath(ld_prune_dir, "shard-$shard", "execution")
    files = readdir(execution_dir)
    fam_file = files[findfirst(endswith(".fam"), files)]
    push!(groups, splitext(splitext(fam_file)[1])[1])
end
@test groups == expected_groups

# Test LOCO PCA
scatter_dirs = filter(x -> occursin("call-Scatter", x), readdir(results_dir, join=true))
loco_pca_dir = findfirst(
    dir_name ->  dir_contains_subdir(dir_name, "call-loco_pca"),
    scatter_dirs
)
loco_pca_dir = scatter_dirs[loco_pca_dir]
## One PCA per (group, chromosome) pair = 5 * 3 = 15
## These are ordered by group and chromosome
pca_groups_and_chrs = Set([])
for group_shard in [0, 1, 2, 3, 4, 5]
    subdir = only(readdir(joinpath(loco_pca_dir, "shard-$group_shard"), join=true))
    subdir = joinpath(only(readdir(subdir, join=true)), "call-loco_pca")
    for chr_shard in [0, 1, 2]
        execution_dir = joinpath(subdir, "shard-$chr_shard", "execution")
        files = readdir(execution_dir)
        eigenvec_file = files[findfirst(endswith("eigenvec"), files)]
        _, ancestry, phenotype, chr, _ = split(eigenvec_file, ".")
        push!(pca_groups_and_chrs, ("$ancestry.$phenotype", chr))
    end
end
@test pca_groups_and_chrs == Set(Iterators.product(expected_groups, ["chr1_out", "chr2_out", "chr3_out"]))

# Test merge covariates and PCs
covariates_and_pcs_dir = joinpath(results_dir, "call-merge_covariates_and_pcs")
merged_covariates_groups = Set([])
for group_shard in 0:5
    execution_dir = joinpath(covariates_and_pcs_dir, "shard-$group_shard", "execution")
    files = readdir(execution_dir)
    merged_covariates_filename = files[findfirst(endswith("merged_covariates_and_pcs.tsv"), files)]
    ancestry, phenotype, _ = split(merged_covariates_filename, ".")
    push!(merged_covariates_groups, "$ancestry.$phenotype")
    covariates_and_pcs = CSV.read(joinpath(execution_dir, merged_covariates_filename), DataFrame)
    for chr in 1:3
        for pc in 1:10
            @test "CHR$(chr)_OUT_PC$(pc)" in names(covariates_and_pcs)
        end
    end
    # Since the dataset is merged with PCs which are computed from individuals with no missing data, there is no missing data for the phenotype of interest
    @test "NA" ∉ covariates_and_pcs[!, phenotype] # missing are coded as NA
    @test eltype(covariates_and_pcs[!, phenotype]) == Int
end
@test merged_covariates_groups == expected_groups

# Test REGENIE Step 1
regenie_step_1_dir = joinpath(results_dir, "call-regenie_step_1")
regenie_groups = Set([])
for group_shard in 0:5
    execution_dir = joinpath(regenie_step_1_dir, "shard-$group_shard", "execution")
    files = readdir(execution_dir)
    pred_list = files[findfirst(endswith(".step1_pred.listrelative"), files)]
    ancestry, phenotype, _ = split(pred_list, ".")
    push!(regenie_groups, "$ancestry.$phenotype")
end
@test regenie_groups == expected_groups

# Test REGENIE Step 2 / finemapping
top_regenie_step_2_dir = findfirst(
    dir_name ->  dir_contains_subdir(dir_name, "call-regenie_step_2"),
    scatter_dirs
)
top_regenie_step_2_dir = scatter_dirs[top_regenie_step_2_dir]
regenie_step_2_groups = Set([])
results_expected_cols = [
        "CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ", "N", "TEST", "BETA", "SE", "CHISQ", "LOG10P", "EXTRA"
    ]
for group_shard in 0:5
    subdir = only(readdir(joinpath(top_regenie_step_2_dir, "shard-$group_shard"), join=true))
    subdir = only(readdir(subdir, join=true))
    regenie_step_2_dir = joinpath(subdir, "call-regenie_step_2")
    for chr_shard in 0:2
        execution_dir = joinpath(regenie_step_2_dir, "shard-$chr_shard", "execution")
        files = readdir(execution_dir)
        step_2_results_file = files[findfirst(endswith(".regenie"), files)]
        ancestry, phenotype, chr_pheno, _ = split(step_2_results_file, ".")
        chr = split(chr_pheno, "_")[1]
        push!(regenie_step_2_groups, ("$ancestry.$phenotype", chr))
        step_2_results = CSV.read(joinpath(execution_dir, step_2_results_file), DataFrame)
        @test names(step_2_results) == results_expected_cols
        @test nrow(step_2_results) > 0
    end
    finemapping_dir = joinpath(subdir, "call-finemapping")
    for chr_shard in 0:2
        execution_dir = joinpath(finemapping_dir, "shard-$chr_shard", "execution")
        fp_files = filter(endswith(".tsv"), readdir(execution_dir))
        @test length(fp_files) == 2
    end
end
@test regenie_step_2_groups == Set(Iterators.product(expected_groups, ["chr1", "chr2", "chr3"]))

# Test Merged gwas results
merge_chr_results_dir = joinpath(results_dir, "call-merge_gwas_group_chr_results")
merged_results_groups = Set([])
for group_shard in 0:5
    execution_dir = joinpath(merge_chr_results_dir, "shard-$group_shard", "execution")
    files = readdir(execution_dir)
    gwas_merged_results_file = files[findfirst(endswith("gwas.tsv"), files)]
    ancestry, phenotype, _ = split(gwas_merged_results_file, ".")
    push!(merged_results_groups, "$ancestry.$phenotype")
end
@test merged_results_groups == expected_groups

# Test Merged finemapping results
merge_chr_results_dir = joinpath(results_dir, "call-merge_fp_group_chr_results")
merged_results_groups = Set([])
for group_shard in 0:5
    execution_dir = joinpath(merge_chr_results_dir, "shard-$group_shard", "execution")
    files = readdir(execution_dir)
    gwas_merged_results_file = files[findfirst(endswith("finemapping.tsv"), files)]
    ancestry, phenotype, _ = split(gwas_merged_results_file, ".")
    push!(merged_results_groups, "$ancestry.$phenotype")
end
@test merged_results_groups == expected_groups

# Test Plots
plots_dir = joinpath(results_dir, "call-gwas_group_plots")
plots_groups = Set([])
for group_shard in 0:5
    execution_dir = joinpath(plots_dir, "shard-$group_shard", "execution")
    files = readdir(execution_dir)
    plot_files = filter(endswith(".png"), files)
    @test length(plot_files) >= 2
    _, _, ancestry, phenotype, _ = split(first(plot_files), ".")
    push!(plots_groups, "$ancestry.$phenotype")
end
@test plots_groups == expected_groups

# Test Meta-analysis
meta_analysis_dir = joinpath(results_dir, "call-meta_analyse", "execution")
meta_analysis_files = filter(endswith(".tsv"), readdir(meta_analysis_dir))
meta_analysed_phenotypes = Set(String[])
for file in meta_analysis_files
    meta_results = CSV.read(joinpath(meta_analysis_dir, file), DataFrame; delim="\t")
    @test length(unique(meta_results.ID)) == length(meta_results.ID)
    @test nrow(meta_results) > 0
    push!(meta_analysed_phenotypes, split(file, ".")[2])
    @test all(meta_results.LOG10P .>= 0.0)
end
@test meta_analysed_phenotypes == Set(["SEVERE_PNEUMONIA", "SEVERE_COVID_19"])

# Test Meta-analysis finemapping
meta_fp_analysis_dir = joinpath(results_dir, "call-merge_fp_meta_chr_results")
for shard in [0, 1] # 2 shards for 2 phenotypes
    execution_dir = joinpath(meta_fp_analysis_dir, "shard-$shard", "execution")
    @test findfirst(endswith("finemapping.tsv"), readdir(execution_dir)) !== nothing
end

# Test Meta-analysis plots
meta_plots_dir = joinpath(results_dir, "call-gwas_meta_plots")
for shard in [0, 1] # 2 shards for 2 phenotypes
    execution_dir = joinpath(meta_plots_dir, "shard-$shard", "execution")
    files = readdir(execution_dir)
    plot_files = filter(endswith(".png"), files)
    @test length(plot_files) == 2
end

end

true