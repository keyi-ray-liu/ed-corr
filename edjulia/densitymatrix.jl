



cdagup(basis_dict:: Dict, ::Electron, td::TwoD, site) = _cdag(basis_dict, Electron(), site, 0, 0 )
cdagdn(basis_dict:: Dict, ::Electron, td::TwoD, site) = _cdag(basis_dict, Electron(), site, 0, td.L )
cdag(args...) = cdagup(args...) + cdagdn(args...)


cup(args...) = cdagup(args...)'
cdn(args...) = cdagdn(args...)'
c(args...) = cdag(args...)'




corr_up( basis_dict, Par, td, site1, site2) =  cdagup( basis_dict, Par, td, site1) * cup(basis_dict, Par, td, site2)
corr_dn( basis_dict, Par, td, site1, site2) = cdagdn( basis_dict, Par, td, site1) * cdn(basis_dict, Par, td, site2) 



function _cdag(basis_dict :: Dict,  ::Electron, site::Int, local_shift::Int, spin_shift::Int  )

    dim = length(keys(basis_dict))
    M = spzeros(dim, dim)
    # this represents cdag at site $site

    for (basis :: Tuple, ind1 :: Int) in basis_dict

        site += local_shift

        if basis[site + spin_shift] == 0
            newbasis =  collect(basis)
            newbasis[site + spin_shift] = 1

            newkey = Tuple(newbasis)

            if newkey ∈ keys(basis_dict)
                ind2 = basis_dict[newkey]
                M[ind2, ind1] = 1
            end 

        end 


    end 

    return M


end 






function lindbladian!(du, u, p, t)

    h = p[1]
    inj = p[2]
    dep = p[3]
    γi = p[4]
    γd = p[5]

    du .= -1im *  commutator(h, u) 
    du .+= γi * ( inj * u * inj' - 1/2 * anticommutator( inj' * inj, u))
    du .+= γd * ( dep * u * dep' - 1/2 * anticommutator( dep' * dep, u))

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
        du += γs[i] * ( op * u * op' - 1/2 * anticommutator( op' * op, u))
    end 

    return du
end


function gen_ρ(qn, Par, Geo)

    basis = gen_basis(qn, Par, Geo)
    dim = length(basis) 
    
    ρ = zeros(dim, dim) .* (1.0 + 0.0im)
    
    ρ[1, 1] = 1

    @show typeof(ρ)
    return ρ
end 


function nup(basis_dict :: Dict, ::Electron, :: Geometry, site)
    return _occ(basis_dict, Electron(), site, 0 )
end 

function ndn(basis_dict :: Dict, ::Electron, Geo:: Geometry, site)
    return _occ(basis_dict, Electron(), site, Geo.L )
end 


function _occ(basis_dict :: Dict,  ::Electron, site::Int,  spin_shift::Number  )

    dim = length(keys(basis_dict))
    M = spzeros(dim, dim)
    # this represents cdag at site $site

    for (basis :: Tuple, ind1 :: Int) in basis_dict
        M[ind1, ind1] = basis[ site + spin_shift] 
    end 

    return M


end 
