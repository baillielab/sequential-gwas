module SequentialGWAS

using CSV
using DataFrames
using DelimitedFiles
using ArgParse
using Random

include("strand_flip.jl")
include("mock.jl")
include("cli.jl")

export julia_main

end
