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
        mkpath( filestr )
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
        mkpath( filestr )
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




function _odesolve(qn::QN, Par :: Fermion, Geo :: Geometry, Coul :: Coulomb, bias :: Bias, ρ0, op :: InjDep;)

    h = gen_ham(qn, Par, Geo, Coul, bias )
    basis = gen_basis(qn, Par, Geo)
    basis_dict = Dict( b => i for (i, b) in enumerate(basis))

    inj_source = cdag(basis_dict,  Geo, op.source_site)
    dep_source = c(basis_dict,  Geo, op.source_site)

    inj_drain= cdag(basis_dict,  Geo, op.drain_site)
    dep_drain = c(basis_dict,  Geo, op.drain_site)

    ops = [inj_source, dep_source, inj_drain, dep_drain ]
    γs = [ op.γ_inj_source, op.γ_dep_source, op.γ_inj_drain, op.γ_dep_drain]

    p = [h, ops, γs]
    tspan = (0, 500)

    
    prob = ODEProblem(lindbladian, ρ0, tspan, p)
    sol = solve(prob, Tsit5(), reltol = 1e-8, abstol = 1e-8)

    return sol

end 



function _odesolve(h, basis_dict, Geo :: Geometry, ρ0, op :: InjDep; start = 0, fin = 500)

    #   ---------- up -------------

    inj_source_up = cdagup(basis_dict, Geo, op.source_site)
    dep_source_up = cup(basis_dict,  Geo, op.source_site)

    inj_drain_up = cdagup(basis_dict, Geo, op.drain_site)
    dep_drain_up = cup(basis_dict,  Geo, op.drain_site)

   # ---------- dn -------------

    inj_source_dn = cdagdn(basis_dict, Geo, op.source_site)
    dep_source_dn = cdn(basis_dict,  Geo, op.source_site)

    inj_drain_dn = cdagdn(basis_dict,  Geo, op.drain_site)
    dep_drain_dn = cdn(basis_dict,  Geo, op.drain_site)


    ops = [
        inj_source_up, dep_source_up, inj_drain_up, dep_drain_up, 
    inj_source_dn, dep_source_dn, inj_drain_dn, dep_drain_dn, 
    ]

    γs = [ 
        op.γ_inj_source, op.γ_dep_source, op.γ_inj_drain, op.γ_dep_drain,
     op.γ_inj_source, op.γ_dep_source, op.γ_inj_drain, op.γ_dep_drain,
     ]

    p = [h, ops, γs]
    tspan = (start, fin)

    
    prob = ODEProblem(lindbladian!, ρ0, tspan, p)

    #method = lsoda()  # recommended for large system but not complex
    #method = DP8()  # supposedly stable memory wise
    #method = VCABM()
    #method = Tsit5()  # out of memory for 500? 
    #method = Vern8()
    method = DP5()

    @time sol = solve(prob, method, reltol = 1e-5, abstol = 1e-5, 
    progress = true, progress_steps = 1
    #saveat=tstep
    )

    #@show sizeof(sol.u) / 2^30

    return sol

end 



function odesolve(qn::QN, Par :: Electron, Geo :: Geometry, Coul :: Coulomb, bias :: Bias, ρ0, op :: InjDep; chunks = 1, start = 0, fin = 500, filestr = "")

    h = gen_ham(qn, Par, Geo, Coul, bias )
    basis = gen_basis(qn, Par, Geo)
    basis_dict = Dict( b => i for (i, b) in enumerate(basis))

    stops = collect(range(start, fin, length = chunks + 1))

    for ind in eachindex(stops[1:end - 1])

        cur_start = stops[ind]
        cur_fin = stops[ind + 1]

        cur_sol = _odesolve(h, basis_dict, Geo, ρ0, op; start = cur_start, fin = cur_fin)

        expectation(Not_conserved(), cur_sol, Par, Geo, Occupation(); filestr = filestr)
        expectation(Not_conserved(), cur_sol, Par, Geo, Current(); filestr = filestr)

        ρ0 = cur_sol.u[end]

        
        
    end 

    @info "memory usage: $(Sys.maxrss()/2^30) GB"

end 
