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

include("resources.jl")
include("one_time_checks.jl")
include("mock.jl")
include("cli.jl")
include("read_write.jl")
include("report.jl")


export julia_main

end
