function dpt_test()

    NW = 1
    S = 4
    U = 5.0
    vs = 3/4
    μ = 1/4
    krylovdim = 20
    maxiter = 100
    tol = 1e-10
    check_exist = false

    dpt = DPT(;NW= NW, S = S, vs = vs)
    initdpt = DPT(;NW= NW, S = S, vs = 0.0)
    systag = get_systag(dpt)

    Coul = Coulomb(U, 0.0, 0.0, 0.0, 0.0)

    γ = 0.9
    μL = [ μ for _ in 1: (NW + div(S, 2))]
    μR = [ -μ for _ in (NW + div(S, 2) + 1): NW * 2 + S]
    DD = [0, 0]
    DDinit = [-1e3, 1e3]

    initbias = vcat([0 for _ in 1: NW * 2 + S], DDinit)
    initbias = Bias(initbias)

    bias = vcat( μL, μR, DD)
    bias = Bias(bias)


    state = Tuple([0 for _ in 1:dpt.L])
    injdep = InjDep(1, NW * 2 + S, γ, 0, 0, γ)
    solver = "exp"
    tfin = 20
    dt = 1/8

    top = "$(homedir())/Desktop/PROJECT_DPT/EDMarkovian/bench/"
    filestr = gen_file(top; 
        U = U,
        sys = systag,
        NW = NW,
        vs = vs,
        tfin = tfin,
        injs = injdep.γ_inj_source,
        deps = injdep.γ_dep_source,
        injd = injdep.γ_inj_drain,
        depd = injdep.γ_dep_drain,
        state = join(state, ""),
        solver = solver,
        dt = dt
    )

    @show state
    @show filestr

    # we solve for the GS
    cf = Fermion(div(dpt.L, 2))
    smallbasis = gen_basis(Conserved(), cf, dpt)
    largebasis = gen_basis(Not_conserved(), Fermion(), dpt)
    w, U = GS(Conserved(), cf, initdpt, Coul, initbias)

    occ = expectation( Conserved(), cf, initdpt, U, Occupation())

    @show occ, sum(occ)

    ρ = gen_ρ(smallbasis, largebasis, U)

    @show @info "memory usage before set up : $(Sys.maxrss()/2^30) GB"

    if !ispath(filestr * "times") || !check_exist
        
        @time times, CC = expsolve(Not_conserved(), Fermion(), dpt, Coul, bias, ρ , injdep; filestr = filestr, start = 0, fin = tfin, dt = dt, krylovdim = krylovdim, maxiter = maxiter, tol = tol, writeateverystep = true)

    else
        @info "data exists! skip cal"
    end 


    plt = plot(times, real.(CC[:, end - 1, end - 1]))
    display(plt)
    return nothing
end 
