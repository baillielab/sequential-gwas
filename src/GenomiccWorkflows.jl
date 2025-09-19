module GenomiccWorkflows

using CSV
using DataFrames
using DelimitedFiles
using ArgParse
using Random
using HTTP
using CodecZlib
using ZipFile
using BufferedStreams
using Downloads
using PackageCompiler
using Tables
using CairoMakie
using Statistics
using Base.Threads
using Literate
using MarkdownTables
using YAML
using GeneticsMakie
using JSON
using MLJBase
using Base.Threads
using Dates
using RCall
using PGENFiles
using MLJTransforms

# Legacy files
include("one_time_checks.jl")
include("mock.jl")
# Common files
include("ancestry.jl")
include("read_write.jl")
# Not yet categorised files
include("resources.jl")
include("cli.jl")
include("report.jl")
include("qc_from_kgp.jl")
include("relatedness.jl")
include("pca_qc.jl")
include("gvcf_genotyping.jl")
# Imputation Files
include("imputation.jl")
# Merging UKB and GenOMICC Files
include("merge_ukb_genomicc.jl")
# GWAS files
include("covariates.jl")
include("gwas.jl")
include("fine_mapping.jl")

export julia_main

end
