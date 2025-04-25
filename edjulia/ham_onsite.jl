function onsite!(basis_dict :: Dict, M :: SparseMatrixCSC, par::Particle, bias :: Bias, Geo :: Geometry)

    for (basis:: Tuple, ind::Int) in basis_dict

        M[ind, ind] += _onsite(basis, par, bias, Geo.L)
    end 

    return M

end 


function _onsite(basis :: Tuple, :: Electron, bias::Bias, L::Int)

    return sum(bias.val .* basis[1:L ] + bias.val .* basis[L + 1:end])

end 


function _onsite(basis :: Tuple, :: Fermion, bias::Bias, L::Int)

    return sum(bias.val .* basis[1:L ])

end 