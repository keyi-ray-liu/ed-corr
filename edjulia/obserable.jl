struct Occupation <: Observable
end 



function expectation(qn::QN, Par::Particle, Geo::Geometry, U, Obs::Observable)
    basis = gen_basis(qn, Par, Geo)

    basis_mat = hcat([collect(arr) for arr in basis]...)
    basis_mat .-= 1

    expectation(basis_mat, Par, Geo, U, Obs)

end 

function expectation(basis_mat :: Union{Matrix, SparseMatrixCSC}, Par::Electron,  Geo :: Geometry , U, ::Occupation)

    L, total = size(basis_mat)

    Up = basis_mat[1:div(L, 2), :]
    Dn = basis_mat[div(L, 2) + 1:end, :]

    occup = []
    occdn = []


    for col in 1:size(U)[2]
        
        v = U[:, col]

        v = abs2.(v)

        upval = Up * v
        dnval = Dn * v

        append!(occup, [upval])
        append!(occdn, [dnval])

    end 

    open( "ref/" * get_name(Geo) * get_name(Par) * "occup", "w") do io
        writedlm(io, occup)
    end

    open( "ref/" * get_name(Geo) * get_name(Par) * "occdn", "w") do io
        writedlm(io, occdn)
    end

    return occup, occdn

end 