struct Occupation <: Observable
end 

struct Current <: Observable end


struct Correlation <: Observable
    hopdict :: Dict
end 

function expectation(qn::QN, Par::Particle, Geo::Geometry, U, ::Occupation; filestr="", save=true)
    basis = gen_basis(qn, Par, Geo)

    basis_mat = hcat([collect(arr) for arr in basis]...)
    basis_mat .-= 1
    occup, occdn = _expectation(basis_mat, Par, Geo, U, Occupation(); filestr=filestr, save=save)

    return occup, occdn
end 

function _expectation(basis_mat :: Matrix, ::Electron,  Geo :: Geometry , V, ::Occupation; filestr="", save=true)

    L = size(basis_mat, 1)

    Up = basis_mat[1:div(L, 2), :]
    Dn = basis_mat[div(L, 2) + 1:end, :]


    V = abs2.(V)
    occup = Up * V
    occdn = Dn * V


    try
        mkdir( filestr )
    catch
    end

    occup, occdn = sanitizer(Occupation(), Geo, occup, occdn)

    if save
        open( "$(filestr)occup", "w") do io
            writedlm(io, round.(occup, sigdigits=5))
        end

        open( "$(filestr)occdn", "w") do io
            writedlm(io, round.(occdn, sigdigits=5))
        end
    end

    #@show sum(occup'[:, 3:end], dims=2)
    return occup, occdn

end 





sanitizer(::Occupation, ::Geometry, occup, occdn) = occup', occdn'
function sanitizer(::Occupation, sd::SD, occup, occdn)

    occup = occup'
    occdn = occdn'

    occup = hcat( occup[ :, 1 :sd.S.L], occup[:, sd.S.L + sd.D.L + 1: end], occup[:, sd.S.L + 1 : sd.S.L + sd.D.L])

    occdn = hcat( occdn[ :, 1 :sd.S.L], occdn[:, sd.S.L + sd.D.L + 1: end], occdn[:, sd.S.L + 1 : sd.S.L + sd.D.L])

    return occup, occdn
end 



function expectation(qn::QN, Par::Electron, Geo::Geometry, U, corr::Correlation; filestr="", save=true)
    basis = gen_basis(qn, Par, Geo)

    basis_dict = Dict( b => i for (i, b) in enumerate(basis))

    
    hopdict = corr.hopdict 

    @show hopdict
    

    corr_mat_up = spzeros( size(basis, 1), size(basis, 1))
    corr_mat_up = _hopping(basis_dict, corr_mat_up, Par, hopdict, 1.0, 0, Geo.L; dn=false, herm=false)


    corr_mat_dn = spzeros( size(basis, 1), size(basis, 1))
    corr_mat_dn = _hopping(basis_dict, corr_mat_dn, Par, hopdict, 1.0, 0, Geo.L; up=false, herm=false)

    corr = _expectation(corr_mat_up, corr_mat_dn, Par, Geo, U, corr; filestr=filestr, save=save)

    return corr


end 

function _expectation(corr_mat_up :: SparseMatrixCSC, corr_mat_dn :: SparseMatrixCSC, ::Electron,  Geo :: Geometry , V, ::Correlation; filestr="", save=true)

    corr = zeros( Complex, size(V, 2), 2)

    for (i, v) in enumerate(eachcol(V))


        corr_up = dot(conj( v'), corr_mat_up, v)
        corr_dn = dot(conj( v'), corr_mat_dn, v)

        corr[i, :] = [corr_up, corr_dn]
    end 


    try
        mkdir( filestr )
    catch
    end


    if save
        open( "$(filestr)corr", "w") do io
            writedlm(io, round.(corr, sigdigits=5))
        end

    end 

    #@show sum(occup'[:, 3:end], dims=2)
    return corr

end 


function expectation(qn::QN, Par::Electron, sd::SD, U, current::Current; filestr="", save=true)


    corr_s = expectation(qn, Par, sd, U, Correlation(sd.S.hopdict); save=false)
    corr_d = expectation(qn, Par, sd, U, Correlation(sd.D.hopdict); save=false)

    corr_s = 2  * imag.(corr_s)
    corr_d = - 2  * imag.(corr_d)

    current = hcat( corr_s[:, 1], corr_d[:, 1], corr_s[:, 2], corr_d[:, 2])

    if save

        open( "$(filestr)currentSD", "w") do io
            writedlm(io, round.(current, sigdigits=5) )
        end

    end 

    return current

end 


