using LinearAlgebra
using Random
using Dates
using DelimitedFiles
using Distributed
using BenchmarkTools
using Combinatorics
using ArnoldiMethod
using SparseArrays

include("structs.jl")
include("obserable.jl")
include("basis.jl")
include("ham_hopping.jl")
include("hamiltonian.jl")
include("tasks.jl")
include("solver.jl")
include("ham_coulomb.jl")
include("ham_hubbard.jl")



task()






# as of Feb 20, work with ED code (Fortran), Fermion, no self NE, 
# work with ED with Electron after JW fix (F_i vs F_j)

# -2.916056032356521
# -4.620045112862984
# -2.8284271247460206