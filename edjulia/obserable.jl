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


function expectation(qn ::QN, sol::ODESolution, Par::Electron, Geo:: Geometry, ::Occupation; filestr="", save=false, if_plot = false)

    Nups, Ndns = _expectation(qn, sol, Par, Geo, Occupation(), if_plot)


    if save

        open( "$(filestr)occup", "w") do io
            writedlm(io, round.(Nups, sigdigits=5) )
        end

        open( "$(filestr)occdn", "w") do io
            writedlm(io, round.(Ndns, sigdigits=5) )
        end
    end 

end 


function expectation(qn ::QN, sol::ODESolution, Par::Electron, Geo:: Geometry, ::Current; filestr="", save=false, if_plot = false)

    cur = _expectation(qn, sol, Par, Geo, Current(), if_plot)


    if save

        open( "$(filestr)current", "w") do io
            writedlm(io, round.(cur, sigdigits=5) )
        end

    end 

end 




function _expectation(qn ::QN, sol::ODESolution, Par::Electron, Geo:: Geometry, ::Occupation, if_plot :: Bool)

    basis = gen_basis(qn, Par, Geo)
    basis_dict = Dict( b => i for (i, b) in enumerate(basis))

    #v = collect(basis)
    #@show v[1], v[2], v[17]

    Nupops = [nup(basis_dict, Electron(), Geo, i) for i in 1:Geo.L]
    Ndnops = [ndn(basis_dict, Electron(), Geo, i) for i in 1:Geo.L]

    Nups = []
    Ndns = []


    for ρ in sol.u

        @show count( !=(0), ρ)
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

    if if_plot
        plt = plot()
        
        for (i, obs) in enumerate(eachcol(Nups))
            plot!(plt, sol.t, obs; label = "n↑$(i)", xlabel = "time", title = "", ylabel = "⟨n̂⟩")
        end 

        # for (i, obs) in enumerate(eachcol(Ndns))
        #     plot!(plt, sol.t, obs; label = "n↓$(i)", xlabel = "time", title = "", ylabel = "⟨n̂⟩")
        # end 

        display(plt)
    end 

    return Nups, Ndns

end 


function _expectation(qn ::QN, sol::ODESolution, Par::Electron, td:: TwoD, ::Current, if_plot :: Bool)

    basis = gen_basis(qn, Par, td)
    basis_dict = Dict( b => i for (i, b) in enumerate(basis))

    L = td.L
    X = td.X
    
    cur_sup = corr_up(basis_dict, Par, td, 1, 2) + corr_up(basis_dict, Par, td, 1, 1 + X)
    cur_sdn = corr_dn(basis_dict, Par, td, 1, 2) + corr_dn(basis_dict, Par, td, 1, 1 + X)

    cur_dup = corr_up(basis_dict, Par, td, L - X, L) + corr_up(basis_dict, Par, td, L - 1, L)
    cur_ddn = corr_dn(basis_dict, Par, td, L - X, L) + corr_dn(basis_dict, Par, td, L - 1, L)


    # cur_sup = cup(basis_dict, Par, td, 1) * cdagup(basis_dict, Par, td, 1) 
    # cur_sdn = cdn(basis_dict, Par, td, 1) * cdagdn(basis_dict, Par, td, 1)
    # cur_dup = cup(basis_dict, Par, td, L) * cdagup(basis_dict, Par, td, L)
    # cur_ddn = cdn(basis_dict, Par, td, L) * cdagdn(basis_dict, Par, td, L)

    curops = [cur_sup, cur_dup, cur_sdn, cur_ddn]
    curs = []
    
    for ρ in sol.u

        #@show ρ
        cur = [ _expectation(ρ, curop) for curop in curops]
        
        append!(curs, [cur])
    end 

    #@show Nups

    curs = vectomat(curs)

    curs = imag.(curs)


    if if_plot
        plt = plot()
        

        for (i, obs) in enumerate(eachcol(curs))
            plot!(plt, sol.t, obs; label = "cur $(i)", xlabel = "time", title = "", ylabel = "I")
        end 

        display(plt)
    end 

    return curs

end 




function _expectation(ρ, op :: SparseMatrixCSC)

    O = tr(ρ * op)
    return O

end 


