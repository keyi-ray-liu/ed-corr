using LinearAlgebra
using Random
using Dates
using DelimitedFiles
using Distributed
using BenchmarkTools
using Combinatorics
using ArnoldiMethod
#using KrylovKit
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
include("ham_onsite.jl")



#task()
#time_evolve()
#GS()
scan()




# as of Feb 20, work with ED code (Fortran), Fermion, no self NE, 
# work with ED with Electron after JW fix (F_i vs F_j)

# -2.916056032356521
# 3Up4DnU100.0 -4.620045112862984
# 4Up4DnU100.0 -2.8284271247460206   delta = -1.79161798
