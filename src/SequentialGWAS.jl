module SequentialGWAS

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

include("resources.jl")
include("one_time_checks.jl")
include("mock.jl")
include("cli.jl")
include("read_write.jl")
include("report.jl")
include("qc_from_kgp.jl")
include("relatedness.jl")
include("pca_qc.jl")
include("ancestry.jl")
include("covariates.jl")

export julia_main

end
