

function gen_ham(qn::QN, Sp :: Fermion, Geo :: Geometry, Coul :: Coulomb)

    basis = gen_basis(qn, Sp, Geo)

    dim = length(basis)
    M = zeros(dim, dim)

    basis_dict = Dict( b => i for (i, b) in enumerate(basis))
    M = hopping!(basis_dict, M, Sp, Geo)
    M = coulomb!(basis_dict, M, Sp, Coul)

    #@show M
    return M
end 