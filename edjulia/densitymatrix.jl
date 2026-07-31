
cdag(basis_dict:: Dict, geo::Geometry, site) = _cdag(basis_dict, site, geo.L)
c(args...) = cdag(args...)'

cdagup(basis_dict:: Dict,  geo::Geometry, site) = _cdag(basis_dict,  site, 0, geo.L)
cdagdn(basis_dict:: Dict,  geo::Geometry, site) = _cdag(basis_dict,  site, geo.L, geo.L)


cup(args...) = cdagup(args...)'
cdn(args...) = cdagdn(args...)'



corr_up( basis_dict,  geo, site1, site2) = cdagup(basis_dict, geo, site1) * cup(basis_dict, geo, site2)
corr_dn( basis_dict,  geo, site1, site2) = cdagdn(basis_dict, geo, site1) * cdn(basis_dict, geo, site2)


corr( basis_dict, geo, site1, site2) = cdag(basis_dict, geo, site1) * c(basis_dict, geo, site2)
# Fermionic cdag
function _cdag(basis_dict :: Dict,  site::Int,  :: Int )

    dim = length(keys(basis_dict))
    M = spzeros(dim, dim)
    # this represents cdag at site $site

    for (basis :: Tuple, ind1 :: Int) in basis_dict

        if basis[site] == 0
            newbasis =  collect(basis)
            newbasis[site] = 1

            newkey = Tuple(newbasis)

            if newkey ∈ keys(basis_dict)
                ind2 = basis_dict[newkey]
                M[ind2, ind1] = jw(Fermion(), basis, site) * 1

            end 

        end 

    end 

    return M


end 


function _cdag(basis_dict :: Dict,  site::Int, spin_shift::Int, L :: Int )

    dim = length(keys(basis_dict))
    M = spzeros(dim, dim)
    # this represents cdag at site $site

    Sp = spin_shift == 0 ? Up() : Dn()

    for (basis :: Tuple, ind1 :: Int) in basis_dict


        if basis[site + spin_shift] == 0
            newbasis =  collect(basis)
            newbasis[site + spin_shift] = 1

            newkey = Tuple(newbasis)

            if newkey ∈ keys(basis_dict)
                ind2 = basis_dict[newkey]
                M[ind2, ind1] = jw(Sp, basis, site, L) * 1

            end 

        end 


    end 

    return M


end 



function lindbladian!(du, u, p, t)

    h = p[1]

    ops = p[2]
    γs = p[3]

    du .= -1im *  commutator(h, u) 

    for (i, op) in enumerate(ops)
        du .+= γs[i] * ( op * u * op' - 1/2 * anticommutator( op' * op, u))
    end 

    #du .= BlockBandedMatrix(du)
    nothing
end


"""Lindbladian in-place operator. p:
p[1] hamiltonain
p[2] = Tuple of operators
p[3] = Tuple of parameters

len(p[2]) and len(p[3]) must match
"""
function lindbladian(u, p, t)

    h = p[1]

    ops = p[2]
    γs = p[3]

    du = -1im *  commutator(h, u) 

    for (i, op) in enumerate(ops)
        #du += γs[i] * ( op * u * op' - 1/2 * anticommutator( op' * op, u))
        du += γs[i] * ( op * u * op' - 1/2 * anticommutator( op' * op, u))
    end 

    return du
end


function gen_ρ(qn, Par, Geo; state = (0, 0, 0, 0, 0, 0, 0, 0))

    basis = gen_basis(qn, Par, Geo)
    dim = length(basis) 
    
    ρ = zeros(ComplexF64, dim, dim) 

    ind = findfirst( x -> x == state, collect(Iterators.flatten((basis,))))

    ρ[ind, ind] = 1

    return ρ
end 


function gen_ρ(smallbasis, largebasis, U)

    dim = length(largebasis)
    
    ρ = zeros(ComplexF64, dim) 
    ref = collect(Iterators.flatten((largebasis,)))

    for (i, val) in enumerate(U)

        state = smallbasis[i]
        ind = findfirst( x -> x == state, ref)

        #@show ind
        ρ[ind] = val

    end
    
    ρ = ρ * ρ'

    return ρ
end 

n(basis_dict :: Dict, Geo :: Geometry, site) = nup(basis_dict, Geo, site)


function nup(basis_dict :: Dict, :: Geometry, site)
    return _n(basis_dict, site, 0 )
end 

function ndn(basis_dict :: Dict,  Geo:: Geometry, site)
    return _n(basis_dict, site, Geo.L )
end 


nupdn(basis_dict :: Dict,  Geo:: Geometry, site) = nup(basis_dict, Geo, site) * ndn(basis_dict, Geo, site)


function _n(basis_dict :: Dict, site::Int,  spin_shift::Number  )

    dim = length(keys(basis_dict))
    M = spzeros(dim, dim)
    # this represents cdag at site $site

    for (basis :: Tuple, ind1 :: Int) in basis_dict
        M[ind1, ind1] = basis[ site + spin_shift] 
    end 

    return M


end 
