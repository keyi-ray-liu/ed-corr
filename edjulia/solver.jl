
# =======================================================================================================================================
# =======================================================================================================================================
# =======================================================================================================================================
# experimental, benchmark liouvillian solver


struct GenericODE end 
struct LinearODE end
struct SteadyState end


# =======================================================================================================================================
# =======================================================================================================================================
# =======================================================================================================================================

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



function _solver(::GenericODE, ρ0, tspan, p)


    prob = ODEProblem(lindbladian!, ρ0, tspan, p)

    @show typeof(p)
    @show typeof(ρ0)
    #@show typeof(inj_source_up)

    #method = lsoda()  # recommended for large system but not complex
    method = DP8()  # supposedly stable memory wise
    #method = VCABM() # very large systemme?
    #method = Tsit5()  # out of memory for 500? 
    #method = Vern7()
    #method = DP5()
    #method = Rodas5P()
    # method = Rodas4P() # stiff for d tol? not vect
    # method = Vern9() #FBDF() # stiff , high tol, vec?
    # method = RKMK4()

    @time sol = solve(prob, method, reltol = 1e-5, abstol = 1e-5, 
    progress = true#, progress_steps = 1
    #saveat=tstep
    )


    # @time solve(prob, KenCarp47(; linsolve = KrylovJL_GMRES());
    # #save_everystep = false
    # )

    #@show sizeof(sol.u) / 2^30
    return sol

end 



function _solver(::LinearODE, ρ0, tspan, p)


    A = MatrixOperator(ones(2, 2), update_func! = lindbladian!)
    prob = ODEProblem(A, ρ0, tspan, p)
    sol = solve(prob, LieRK4(), dt = 1 / 4)

    return sol

end 

function _solver(::SteadyState, ρ0, tspan, p)
    prob = ODEProblem(lindbladian!, ρ0, tspan, p)

    SSprob = SteadyStateProblem(prob)

    @time sol = solve(SSprob, SSRootfind())

end 


function _gen_ops(h, basis_dict, Geo :: Geometry, op :: InjDep)

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


    # ops = [Matrix(op) for op in ops]
    # h = Matrix(h)


    p = [h, ops, γs]

    return p

end 

function _odesolve(h, basis_dict, Geo :: Geometry, ρ0, op :: InjDep; start = 0, fin = 500, backendstr = "generic")

    #   ---------- up -------------

    p = _gen_ops(h, basis_dict, Geo, op)
    tspan = (start, fin)
    
    if backendstr == "generic"
        backend = GenericODE()

    elseif backendstr == "linear"
        backend = LinearODE()

    elseif backendstr == "steadystate"
        backend = SteadyState()

    else
        error("unknown backend")
    end 

    
    # sparsity detector
    #check_sparse(lindbladian!, ρ0, p)
    sol = _solver(backend, ρ0, tspan, p)

    return sol

end 



function odesolve(qn::QN, Par :: Electron, Geo :: Geometry, Coul :: Coulomb, bias :: Bias, ρ0, op :: InjDep; chunks = 1, start = 0, fin = 500, filestr = "", backendstr = "generic")

    h = gen_ham(qn, Par, Geo, Coul, bias )
    basis = gen_basis(qn, Par, Geo)
    basis_dict = Dict( b => i for (i, b) in enumerate(basis))

    stops = collect(range(start, fin, length = chunks + 1))
    

    for ind in eachindex(stops[1:end - 1])

        cur_start = stops[ind]
        cur_fin = stops[ind + 1]

        cur_sol = _odesolve(h, basis_dict, Geo, ρ0, op; start = cur_start, fin = cur_fin, backendstr = backendstr)

        #expectation(Not_conserved(), cur_sol, Par, Geo, Occupation(); filestr = filestr)
        #expectation(Not_conserved(), cur_sol, Par, Geo, Current(); filestr = filestr)
        expectation(Not_conserved(), cur_sol, Par, Geo, Correlation(); filestr = filestr)

        ρ0 = cur_sol.u[end]

        
        
    end 

    @info "memory usage: $(Sys.maxrss()/2^30) GB"

end 




function expsolve(qn::QN, Par :: Electron, Geo :: Geometry, Coul :: Coulomb, bias :: Bias, ρ0, op :: InjDep; start = 0, fin = 50, dt = 1/4, filestr = "")

    h = gen_ham(qn, Par, Geo, Coul, bias )
    basis = gen_basis(qn, Par, Geo)
    basis_dict = Dict( b => i for (i, b) in enumerate(basis))

    p = _gen_ops(h, basis_dict, Geo, op)
    N = size(basis, 1)

    L = Geo.L
    T = length(start:dt:fin)

    upops = Dict((from, to) => corr_up(basis_dict, Geo, from, to) for from in 1:L for to in from:L)
    dnops = Dict((from, to) => corr_dn(basis_dict, Geo, from, to) for from in 1:L for to in from:L)

    CCups = zeros(ComplexF64, T, L, L)
    CCdns = zeros(ComplexF64, T, L, L)
    
    function L_action(u)

        u = reshape(u, N, N)

        h = p[1]
        ops = p[2]
        γs = p[3]

        du = -1im *  commutator(h, u) 

        for (i, op) in enumerate(ops)
            du .+= γs[i] * ( op * u * op' - 1/2 * anticommutator( op' * op, u))
        end 

        #du .= BlockBandedMatrix(du)
        du = vec(du)

        return du
    end

    #L_map = LinearMap(L_action!, N^2, ismutating=true, ishermitian = false, issymmetric = false)
    ρ = ρ0
    # expectation at 0
    for i in 1:L
        for j in i:L
            upval = _expectation(ρ, upops[(i, j)])
            dnval = _expectation(ρ, dnops[(i, j)])
            CCups[1, i, j] = upval
            CCups[1, j, i] = conj(upval)

            CCdns[1, i, j] = dnval
            CCdns[1, j, i] = conj(dnval)
        end 
    end 

    for (tt, t) in enumerate(start:dt:fin)

        @show tt, t

        #ρ_vec = expv(t, L_map, vec(ρ0))
        @time ρ_vec, info = exponentiate(L_action, dt, vec(ρ))
        @show info

        ρ = reshape(ρ_vec, N, N)
        
        @show size(ρ)
        
        # expectation
        for i in 1:L
            for j in i:L
                upval = _expectation(ρ, upops[(i, j)])
                dnval = _expectation(ρ, dnops[(i, j)])
                CCups[tt + 1, i, j] = upval
                CCups[tt + 1, j, i] = conj(upval)

                CCdns[tt + 1, i, j] = dnval
                CCdns[tt + 1, j, i] = conj(dnval)
            end 
        end 
  
    end 


    T, Nx, Ny = size(CCups)
    try
        mkpath( filestr )
    catch
    end


    open( "$(filestr)time", "w") do io
        writedlm(io, round.(collect(start:dt:fin), sigdigits=5) )
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

    @info "memory usage: $(Sys.maxrss()/2^30) GB"

end 
