using LinearAlgebra
using Random
using Dates
using DelimitedFiles
using Distributed
using BenchmarkTools
using Combinatorics

include("structs.jl")
include("basis.jl")
include("ham_hopping.jl")
include("hamiltonian.jl")
include("tasks.jl")
include("solver.jl")
include("ham_coulomb.jl")

task()

