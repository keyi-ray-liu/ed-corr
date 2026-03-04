struct Occupation <: Observable
end 

struct Current <: Observable end


struct Correlation <: Observable
    hopdict :: Dict
    function Correlation(; hopdict = Dict())
        new(hopdict)
    end
    Correlation(hopdict) = Correlation(; hopdict =hopdict)
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
        mkpath( filestr )
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
# function sanitizer(::Occupation, sd::SD, occup, occdn)

#     occup = occup'
#     occdn = occdn'

#     occup = hcat( occup[ :, 1:sd.S.L], occup[:, (sd.S.L + sd.D.L + 1):end], occup[:, (sd.S.L + 1):(sd.S.L + sd.D.L)])
#     occdn = hcat( occdn[ :, 1:sd.S.L], occdn[:, (sd.S.L + sd.D.L + 1):end], occdn[:, (sd.S.L + 1):(sd.S.L + sd.D.L)])

#     return occup, occdn
# end 



function expectation(qn::QN, Par::Electron, Geo::Geometry, U, corr::Correlation; filestr="", save=true)
    basis = gen_basis(qn, Par, Geo)
    basis_dict = Dict( b => i for (i, b) in enumerate(basis))

    
    hopdict = corr.hopdict 

    @show hopdict
    

    corr_mat_up = spzeros( size(basis, 1), size(basis, 1))
    corr_mat_up = _hopping(basis_dict, corr_mat_up, Par, hopdict, 0; dn=false, herm=false)


    corr_mat_dn = spzeros( size(basis, 1), size(basis, 1))
    corr_mat_dn = _hopping(basis_dict, corr_mat_dn, Par, hopdict, Geo.L; up=false, herm=false)

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
        mkpath( filestr )
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


function expectation(qn::QN, Par::Electron, sd::SD, U, ::Current; filestr="", save=true)


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



function expectation(qn ::QN, sol::ODESolution, Par::Fermion, Geo:: Geometry, ::Occupation; filestr="")

    try
        mkpath( filestr )
    catch
    end

    N = _expectation(qn, sol, Par, Geo, Occupation())

    open( "$(filestr)occ", "a") do io
        writedlm(io, round.(N, sigdigits=5) )
    end

end 


function expectation(qn ::QN, sol::ODESolution, Par::Electron, Geo:: Geometry, ::Occupation; filestr="")

    try
        mkpath( filestr )
    catch
    end

    Nups, Ndns = _expectation(qn, sol, Par, Geo, Occupation())

    open( "$(filestr)occup", "a") do io
        writedlm(io, round.(Nups, sigdigits=5) )
    end

    open( "$(filestr)occdn", "a") do io
        writedlm(io, round.(Ndns, sigdigits=5) )
    end


end 


function expectation(qn ::QN, sol::ODESolution, Par::Particle, Geo:: Geometry, ::Current; filestr="")

    try
        mkpath( filestr )
    catch
    end

    cur = _expectation(qn, sol, Par, Geo, Current())
    
    

    open( "$(filestr)time", "a") do io
        writedlm(io, round.(sol.t, sigdigits=5) )
    end

    open( "$(filestr)current", "a") do io
        writedlm(io, round.(cur, sigdigits=5) )
    end


end 

function expectation(qn ::QN, sol::ODESolution, Par::Particle, Geo:: Geometry, ::Correlation; filestr="")

    try
        mkpath( filestr )
    catch
    end

    CCups, CCdns = _expectation(qn, sol, Par, Geo, Correlation())
    T, Nx, Ny = size(CCups)

    open( "$(filestr)time", "a") do io
        writedlm(io, round.(sol.t, sigdigits=5) )
    end

    h5open("$(filestr)CC.h5", "w") do f

        dRu = create_dataset(f, "REup", Float64, dataspace(T, Nx, Ny))
        dRd = create_dataset(f, "REdn", Float64, dataspace(T, Nx, Ny))
        dIu = create_dataset(f, "IMup", Float64, dataspace(T, Nx, Ny))
        dId = create_dataset(f, "IMdn", Float64, dataspace(T, Nx, Ny))

        write(dRu, real.(CCups))
        write(dRd, real.(CCdns))
        write(dIu, imag.(CCups))
        write(dId, imag.(CCdns))

    end 


end 


function _expectation(qn ::QN, sol::ODESolution, Par::Fermion, Geo:: Geometry, ::Occupation)

    basis = gen_basis(qn, Par, Geo)
    basis_dict = Dict( b => i for (i, b) in enumerate(basis))

    #v = collect(basis)
    #@show v[1], v[2], v[17]

    Nops= [n(basis_dict, Geo, i) for i in 1:Geo.L]
    Ns = []


    for ρ in sol.u

        #@show count( !=(0), ρ)
        #@show ρ
        N = [ _expectation(ρ, Nop) for Nop in Nops]
        
        append!(Ns, [N])
    end 
    
    Ns = vectomat(Ns)
    Ns = real.(Ns)


    return Ns

end 



function _expectation(qn ::QN, sol::ODESolution, Par::Electron, Geo:: Geometry, ::Occupation)

    basis = gen_basis(qn, Par, Geo)
    basis_dict = Dict( b => i for (i, b) in enumerate(basis))

    #v = collect(basis)
    #@show v[1], v[2], v[17]

    Nupops = [nup(basis_dict, Geo, i) for i in 1:Geo.L]
    Ndnops = [ndn(basis_dict, Geo, i) for i in 1:Geo.L]

    Nups = []
    Ndns = []


    for ρ in sol.u

        #@show count( !=(0), ρ)
        #@show ρ
        Nup = [ _expectation(ρ, Nupop) for Nupop in Nupops]
        Ndn = [ _expectation(ρ, Ndnop) for Ndnop in Ndnops]
        
        append!(Nups, [Nup])
        append!(Ndns, [Ndn])
    end 

    #@show Nups

    Nups = vectomat(Nups)
    Ndns = vectomat(Ndns)

    Nups = real.(Nups)
    Ndns = real.(Ndns)

    return Nups, Ndns

end 



source(sd :: SD) = [(sd.scoup, 1, 3)]
drain(sd :: SD) = [(sd.dcoup, 2, sd.L)]


source(td :: TwoD) = [(t, 1, to) for (t, to) in td.hopdict[1]]
drain(td :: TwoD) = [(td.hopdict[td.L - td.X][1][1], td.L - td.X, td.L), (td.hopdict[td.L - 1][1][1], td.L - 1, td.L)]


function _expectation(qn ::QN, sol::ODESolution, Par::Fermion, Geo:: Geometry, ::Current)

    basis = gen_basis(qn, Par, td)
    basis_dict = Dict( b => i for (i, b) in enumerate(basis))

    sourceinds = source(Geo)
    draininds = drain(Geo)

    cur_s = sum([coup * corr(basis_dict, td, from, to) for (coup, from, to) in sourceinds ])
    cur_d = sum([coup * corr(basis_dict, td, from, to) for (coup, from, to) in draininds ])



    curops = [cur_s, cur_d]
    curs = []
    
    for ρ in sol.u

        #@show ρ
        cur = [ _expectation(ρ, curop) for curop in curops]
        append!(curs, [cur])
    end 


    curs = vectomat(curs)
    curs = -2 * imag.(curs)

    return curs

end 


function _expectation(qn ::QN, sol::ODESolution, Par::Electron, Geo:: Geometry, ::Current)

    basis = gen_basis(qn, Par, Geo)
    basis_dict = Dict( b => i for (i, b) in enumerate(basis))

    sourceinds = source(Geo)
    draininds = drain(Geo)

    CCups, CCdns = _expectation(qn, sol, Par, Geo, Correlation())
    curs = []
    
    for ρ in sol.u

        #@show ρ
        cur = [ _expectation(ρ, curop) for curop in curops]
        
        append!(curs, [cur])
    end 

    #@show Nups

    curs = vectomat(curs)
    curs = -2 * imag.(curs)

    return curs

end 


function _expectation(qn ::QN, sol::ODESolution, Par::Electron, Geo:: Geometry, ::Correlation)

    basis = gen_basis(qn, Par, Geo)
    basis_dict = Dict( b => i for (i, b) in enumerate(basis))
    L = Geo.L


    upops = Dict((from, to) => corr_up(basis_dict, Geo, from, to) for from in 1:L for to in from:L)
    dnops = Dict((from, to) => corr_dn(basis_dict, Geo, from, to) for from in 1:L for to in from:L)

    T = length(sol.t)

    CCups = zeros(ComplexF64, T, L, L)
    CCdns = zeros(ComplexF64, T, L, L)

    @show size(CCups)
    
    for (tt, ρ) in enumerate(sol.u)

        for i in 1:L
            for j in i:L
                upval = _expectation(ρ, upops[(i, j)])
                dnval = _expectation(ρ, dnops[(i, j)])
                CCups[tt, i, j] = upval
                CCups[tt, j, i] = conj(upval)

                CCdns[tt, i, j] = dnval
                CCdns[tt, j, i] = conj(dnval)
            end 
        end 

        #@show CCups[tt, :, :]
    end 


    return CCups, CCdns

end 




function _expectation(ρ, op :: SparseMatrixCSC)

    O = tr(ρ * op)
    return O

end 
