
using Pkg

for strs in [
    # "Combinatorics",
    # "ArnoldiMethod",
    # "DifferentialEquations",
    # "SparseArrays",
    # "TerminalLoggers",
    # "LSODA",
    # "ForwardDiff",
    # "BenchmarkTools",
    # "BlockBandedMatrices"
    # "Plots",
    # "SparseConnectivityTracer",
    # "ADTypes",
    # "HDF5"
    # "FileIO"
    # "OrdinaryDiffEqLinear",
    # #"LinearMaps",
    # #"ExponentialUtilities",
    # "KrylovKit"
    ]

    Pkg.add(strs)
end 


#using LinearMaps
using TerminalLoggers#, ProgressLogging
using LinearAlgebra
using Random
using Dates
using DelimitedFiles
using DifferentialEquations
using ForwardDiff
using Distributed
using BenchmarkTools
using Combinatorics: combinations
using ArnoldiMethod
using BlockBandedMatrices
using KrylovKit
using SparseArrays
using LSODA
using HDF5
using Plots
using FileIO
using OrdinaryDiffEqLinear
#using ExponentialUtilities
import SparseConnectivityTracer, ADTypes

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
include("densitymatrix.jl")
include("utils.jl")
include("plotter.jl")
include("scan.jl")
include("test.jl")

ForwardDiff.can_dual(::Type{ComplexF64}) = true

if ARGS == []
#task()
    #time_evolve_test()
    
    #markov_test()
    #gamma_scan(; Geo = FRUSTRATION_4())
    #gamma_scan(; Geo = FIVE())
    gamma_scan()
#GS()
#scan()
    #GG()
    #comb_test()
    #check_ham_construction()
    #test()
    #markov_test()
    #setup_test()

#test()
#gamma_scan()
else

    if ARGS[1] == "GG"
        
        GG()
        #GG(float(ARGS[2]), float(ARGS[3]))

    elseif ARGS[1] == "Gamma"

        gamma_scan()
    end 


end 

#markov_test()
#check_ham_construction()




# as of Feb 20, work with ED code (Fortran), Fermion, no self NE, 
# work with ED with Electron after JW fix (F_i vs F_j)

# -2.916056032356521
# 3Up4DnU100.0 -4.620045112862984
# 4Up4DnU100.0 -2.8284271247460206   delta = -1.79161798


# work without Coulomb
# Coulomb is fixed by removing the factor of 2 on NE
# 

# May 28, with only up, dn, works 
# May 29, try to debug with jordanwiger string. 
# May 30 SPinless dynamics is currently broken



# Fermion 3x3 4 -16.432764915936506
# Electron 

# Jun 12 revamp solver to break into smaller 'chunks', Fermion broken