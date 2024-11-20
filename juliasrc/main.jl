using LinearAlgebra
using Random
using Dates
using DelimitedFiles
using Distributed
using BenchmarkTools

include("systems.jl")
include("solver.jl")
include("tasks.jl")

#SDscan()
#single()
SmallScan()
