using Pkg

for strs in [
    "Plots",
    "Interpolations",
    "BenchmarkTools"
    ]

    Pkg.add(strs)
end 




using LinearAlgebra
using Random
using Dates
using DelimitedFiles
using Distributed
using BenchmarkTools
using Interpolations
#using Plots
#using GLMakie
using Plots

include("utils.jl")
include("systems.jl")
include("solver.jl")
include("tasks.jl")


#SDscan()
#single()
#SmallScan()
#Scalingtest()
#CSScan()
#LengthScale()
manual()