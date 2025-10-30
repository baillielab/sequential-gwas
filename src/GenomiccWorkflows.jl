module GenomiccWorkflows

const ENSEMBL_SERVER = "https://rest.ensembl.org"

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
using Tables
using CairoMakie
using Colors
using Statistics
using Base.Threads
using Literate
using MarkdownTables
using YAML
using JSON
using Base.Threads
using GenomiccUtils


include("read_write.jl")
include("resources.jl")
include("cli.jl")
include("report.jl")
include("qc_from_kgp.jl")
include("relatedness.jl")
include("pca_qc.jl")
include("gvcf_genotyping.jl")
include("covariates.jl")

export julia_main

end
