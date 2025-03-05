

function gen_ham(qn::QN, Par :: Fermion, Geo :: Geometry, Coul :: Coulomb )

    basis = gen_basis(qn, Par, Geo)

    dim = length(basis)
    M = zeros(dim, dim)

    basis_dict = Dict( b => i for (i, b) in enumerate(basis))
    M = hopping!(basis_dict, M, Par, Geo)
    M = coulomb!(basis_dict, M, Par, Coul, Geo)

    #@show M
    return M
end 


function gen_ham(qn::QN, Par :: Electron, Geo :: Geometry, Coul :: Coulomb, bias :: Bias)

    basis = gen_basis(qn, Par, Geo)

    dim = length(basis)
    M = spzeros(dim, dim)
    #M = zeros(dim, dim)

    basis_dict = Dict( b => i for (i, b) in enumerate(basis))


    M = hopping!(basis_dict, M, Par, Geo)
    M = coulomb!(basis_dict, M, Par, Coul, Geo)
    M = hubbard(basis_dict, M, Par, Geo)
    M = onsite!(basis_dict, M, Par, bias, Geo)


    @info "setup complete"
    #@show M
    return M
end 