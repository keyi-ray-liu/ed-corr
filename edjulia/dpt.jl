function dpt_single(; NW = 1, S = 4, U = 2.0, vs = 3/4, μ = 1/4, tfin = 140, tswitch = tfin/2, γ = 0.5, dt = 1/8, top = "$(homedir())/Desktop/PROJECT_DPT/EDMarkovian/bench/" )

    krylovdim = 20
    maxiter = 100
    tol = 1e-10
    check_exist = true

    mid = NW + div(S, 2)

    dpt = DPT(;NW= NW, S = S, vs = vs)
    initdpt = DPT(;NW= NW, S = S, vs = 0.0)
    systag = get_systag(dpt)

    Coul = Coulomb(U, 0.0, 0.0, 0.0, 0.0)

    μL = [ μ for _ in 1: mid]
    μR = [ -μ for _ in mid + 1: NW * 2 + S]
    DD = [0, 0]
    DDinit = [-1e3, 1e3]

    initbias = vcat([0 for _ in 1: NW * 2 + S], DDinit)
    initbias = Bias(initbias)

    bias = vcat( μL, μR, DD)
    bias = Bias(bias)

    @info "initbias: $initbias"
    @info "bias: $bias"

    #state = Tuple([0 for _ in 1:dpt.L])
    injdep = InjDep(1, NW * 2 + S, γ, 0.0, 0.0, γ)
    solver = "exp"
    
    filestr = gen_file(top; 
        U = U,
        sys = systag,
        NW = NW,
        vs = vs,
        tfin = tfin,
        tswitch = tswitch,
        injs = injdep.γ_inj_source,
        deps = injdep.γ_dep_source,
        injd = injdep.γ_inj_drain,
        depd = injdep.γ_dep_drain,
        #state = join(state, ""),
        gamma = γ,
        mu = μ,
        S = S,
        solver = solver,
        dt = dt
    )

    #@show state
    @info "filestr: $filestr"

    if !ispath(filestr * "DONE") || !check_exist

        # if isfile("$(filestr)times")
        #     time = readdlm( "$(filestr)times", Float64)
        #     if time[end] == tfin
        #         @info "data exists! skip cal"
        #         writedlm( "$(filestr)DONE", ["DONE"])
        #         return nothing
        #     end
        # end 


        if isfile( "$(filestr)rho.h5")
            @info "ρ exists, load it"
            ρ = h5open("$(filestr)rho.h5", "r") do file
                read(file["rho"])
            end



            rhotime = readdlm( "$(filestr)rhotime", Float64)
            start = rhotime[end]

            CCtag = start < tswitch ? "S2" : "S3"

            restart = h5open("$(filestr)CC$(CCtag).h5", "r") do file
                read(file["restart"])
            end

            restart += 1
        else
            @info "ρ does not exist, generate it"

            # we solve for the GS
            cf = Fermion(div(dpt.L, 2))
            smallbasis = gen_basis(Conserved(), cf, dpt)
            largebasis = gen_basis(Not_conserved(), Fermion(), dpt)
            w, U = GS(Conserved(), cf, initdpt, Coul, initbias)

            occ = expectation( Conserved(), cf, initdpt, U, Occupation())

            @info "occupation numbers: $occ, sum: $(sum(occ))"
            @info "memory usage before set up : $(Sys.maxrss()/2^30) GB"

            ρ = gen_ρ(smallbasis, largebasis, U)
            start = 0.0
            restart = 0
        end

        CCinds = [ (dpt.L - 1, dpt.L - 1), (mid, mid + 1)]
        CCtags = ["n1", "cur"]
        CCfuns = [ val -> real.(val), val -> -2. * imag.(val)]

        if start < tswitch
            @time ρ = expsolve(Not_conserved(), Fermion(), initdpt, Coul, bias, ρ , injdep; filestr = filestr, start = start, fin = tswitch, dt = dt, krylovdim = krylovdim, maxiter = maxiter, tol = tol, writeateverystep = false, tag = "S2", CCinds = CCinds, CCtags = CCtags, CCfuns = CCfuns, restart = restart)
        end
        
        @time ρ = expsolve(Not_conserved(), Fermion(), dpt, Coul, bias, ρ , injdep; filestr = filestr, start = max(start, tswitch), fin = tfin, dt = dt, krylovdim = krylovdim, maxiter = maxiter, tol = tol, writeateverystep = false,  tag = "S3", CCinds = CCinds, CCtags = CCtags, CCfuns = CCfuns, restart = restart)


        # if ifplot
        #     plt = plot(times, n1, xlabel = "time", ylabel = "⟨n₁⟩", label = "⟨n₁⟩", lw = 2, color =:blue, legend = :topright)
        #     right = twinx(plt)
        #     plot!(right, times, cur,  ylabel = "I", label = "I", lw = 2, color =:red, legend = :bottomright)
        #     display(plt)
        # end 

        open( "$(filestr)DONE", "w") do io
            writedlm(io, ["DONE"])
        end

    else
        @info "data exists! skip cal"
    end 



    return nothing
end 


function dpt_test()

    NWs = [2, 3]
    Us = [2.0, 3.0, 3.5, 4.0, 5.0, 6.0]
    vss = [1/32, 1/8, 3/4]
    γs = [1/10, 1/4, 1/3, 1/2, 5/6]

    for NW in NWs
        Threads.@threads for (U, vs, γ) in collect(Base.product( Us, vss, γs))

            @show (NW, U, vs, γ)

            top = joinpath(pwd(), "DPT", "Jul31INIT")
            dpt_single(; NW = NW, tfin = 256, U = U, γ = γ, vs = vs, top = top)
            
        end 
    end 


end 