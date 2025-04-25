# function solve(h, Geo :: Geometry)

#     @assert h == transpose(h)
#     h = Hermitian(h)
#     eigen_h = @time eigen(h, 1:min(1, length(h)))

#     w = eigen_h.values
#     U = eigen_h.vectors

#     @show w
#     @show length(w)
#     @show w[1], w[end]

#     open( "ref/energy" * get_name(Geo), "w") do io
#         writedlm(io, w)
#     end

# end 



function _solve(h, Geo :: Geometry, Par :: Particle; method ::DiagMethod = full_diag(), nev=min(2, length(h)))

    @show nev
    @assert h == transpose(h)
    
    #eigen_h = @time eigen(h, 1:min(1, length(h)))

    w, U = _diag(h, method, nev)

    @show size(U)

    # w, U, info = eigsolve(h,  nev, :SR)

    # @show info
    # @show size(w), size(U)


    filestr = "ref/$(get_name(Geo))$(get_name(Par))/"

    try
        mkdir( filestr )
    catch
    end

    open( "$(filestr)energy" , "w") do io
        writedlm(io, w)
    end

    return w, U

end 


function _time_evolve(qn::QN, Geo_init::Geometry, Geo_dyna::Geometry, Par::Particle, Coul::Coulomb, bias_dyna::Bias, bias_init::Bias, timecontrol :: TimeControl, obs :: Array{T} where T <: Observable, filestr::String)

    @show Geo_init
    @show Geo_dyna


    times = collect(0:timecontrol.dt : timecontrol.fin)

    h_dyna = gen_ham(qn, Par, Geo_dyna, Coul, bias_dyna )
    h_init = gen_ham(qn, Par, Geo_init, Coul, bias_init)

    method =  size(h_dyna, 1) < 5000 ? full_diag() : arnoldi()

    _, v = _solve(h_init, Geo_init, Par; method = arnoldi(), nev=1)
    w, U = _solve(h_dyna, Geo_dyna, Par; method = method, nev=min(1000, size(v, 1) - 2))

    # we find the coefficients of the all the overlaps, then explicitly time_evolve

    coeffs = vec(conj(transpose(v[:, 1])) * U)
    E = exp.(im * w * times')

    @show size(coeffs)
    @show size(E)

    if sum( abs2, coeffs) < 0.99
        @warn "sum smaller than 0.99"
    end 


    newcoeffs = coeffs' .* E'
    @show size(newcoeffs)

    V = U * newcoeffs'
    @show size(V)


    for ob in obs
        @time expectation(qn, Par, Geo_dyna, V, ob; filestr =filestr)
    end 


    try
        mkdir( filestr )
    catch
    end

    open( "$(filestr)times", "w") do io
        writedlm(io, times)
    end


end 



function _diag(h, ::full_diag, nev::Int)

    k = Matrix(h)
    @time w, U = eigen(k)

    return w[1:nev], U[:, 1:nev]
end 

function _diag(h, ::arnoldi, nev)

    h = Hermitian(h)

    @time eigen_h, history = partialschur(h, nev =nev, tol=1e-7, which=:SR, restarts=1000)

    @show history

    w = real(eigen_h.eigenvalues)
    U = eigen_h.Q

    return w, U
end 


function _odesolve(qn::QN, Par :: Electron, Geo :: Geometry, Coul :: Coulomb, bias :: Bias, ρ0, op :: InjDep)

    h = gen_ham(qn, Par, Geo, Coul, bias )
    basis = gen_basis(qn, Par, Geo)
    basis_dict = Dict( b => i for (i, b) in enumerate(basis))



    inj_source = cdag(basis_dict, Par, Geo, op.source_site)
    dep_source = c(basis_dict, Par, Geo, op.source_site)

    inj_drain= cdag(basis_dict, Par, Geo, op.drain_site)
    dep_drain = c(basis_dict, Par, Geo, op.drain_site)

    ops = [inj_source, dep_source, inj_drain, dep_drain]
    γs = [ op.γ_inj_source, op.γ_dep_source, op.γ_inj_drain, op.γ_dep_drain]

    p = [h, ops, γs]
    tspan = (0, 100)

    
    prob = ODEProblem(lindbladian, ρ0, tspan, p)
    sol = solve(prob, Tsit5(), reltol = 1e-8, abstol = 1e-8)

    return sol

end 