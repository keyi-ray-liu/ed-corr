
function GG_scan(; Geo = nothing)


    
    γ = 0.02
    U = 4.0
    tfin = 500
    dt = 500

    krylovdim = 90
    maxiter = 500
    tol = 1e-11


    U = trunc(U; sigdigits = 5)

    if isnothing(Geo)
        Geo = TwoD(2, 2)
    end 
    #Geo = SD(2, 2; scoup = -0.02, dcoup = -0.02)
    #Par = Fermion() 

    γ = trunc(γ; sigdigits = 5)

    Coul = RefCoul
    bias = Bias([0 for _ in 1:Geo.L])

    state = Tuple([0 for _ in 1:Geo.L * 2])
    injdep = InjDep(1, 4, γ, 0.0 , 0.0, γ)

    top = "/home/keyi-liu/Desktop/Code/Markovian/Mar23GGScan/"
    one_run(; γ = γ, U = U, tfin = tfin, dt = dt,  krylovdim = krylovdim, maxiter = maxiter, tol = tol, Geo = Geo, solver = "exp", Coul = Coul, bias = bias, state = state, injdep = injdep, top = top, plot = false
    )

    return nothing
end 



function gamma_scan(; Geo = nothing)


    
    γs = 10.0 .^ (-2:0.1:3.0) #[0.1, 1.0, 10.0, 100.0]
    Us = [0.0] #[1.0, 10.0, 100.0]
    tfin = 750
    dt = 250

    krylovdim = 90
    maxiter = 500
    tol = 1e-11

    for U in Us

        U = trunc(U; sigdigits = 5)

        if isnothing(Geo)
            Geo = TwoD(2, 2)
        end 
        #Geo = SD(2, 2; scoup = -0.02, dcoup = -0.02)
        #Par = Fermion() 
        Threads.@threads for γ in γs

            γ = trunc(γ; sigdigits = 5)

            Coul = ZeroCoul
            bias = Bias([0 for _ in 1:Geo.L])

            state = Tuple([0 for _ in 1:Geo.L * 2])
            injdep = InjDep(1, 4, γ, 0.0 , 0.0, γ)

            top = "/home/keyi-liu/Desktop/Code/Markovian/Mar13EXP/"

            one_run(; γ = γ, U = U, tfin = tfin, dt = dt,  krylovdim = krylovdim, maxiter = maxiter, tol = tol, Geo = Geo, solver = "exp", Coul = Coul, bias = bias, state = state, injdep = injdep, top = top, plot = false)
        end 
    end 

    return nothing
end 




function one_run(; γ = nothing, U = nothing, tfin = nothing, dt = nothing,  krylovdim = 90, maxiter = 500, tol = 1e-11, Geo = nothing, solver = "exp", Coul = nothing, bias = nothing, state = nothing, injdep = nothing, top = nothing, plot = false)


    U = trunc(U; sigdigits = 5)

    if isnothing(Geo)
        Geo = TwoD(2, 2)
    end 
    #Geo = SD(2, 2; scoup = -0.02, dcoup = -0.02)
    systag = get_systag(Geo)
    #Par = Fermion() 

    γ = trunc(γ; sigdigits = 5)
    Par = Electron(; U = U)

    filestr = gen_file(top; 
        U = Par.U,
        Coul = Coul.ee,
        sys = systag,
        tfin = tfin,
        injs = injdep.γ_inj_source,
        deps = injdep.γ_dep_source,
        injd = injdep.γ_inj_drain,
        depd = injdep.γ_dep_drain,
        GOne = G1,
        GTwo = G2,
        devicebias = devicebias,
        state = join(state, ""),
        solver = solver,
        dt = dt
    )

    if !ispath(filestr * "time")
        ρ = gen_ρ(Not_conserved(), Par, Geo; state = state)

        if solver == "ODE"
            @time odesolve(Not_conserved(), Par, Geo, Coul, bias, ρ , injdep; filestr = filestr, start = 0, fin = tfin, chunks = 1, backendstr = "generic")
        else
            @time expsolve(Not_conserved(), Par, Geo, Coul, bias, ρ , injdep; filestr = filestr, start = 0, fin = tfin, dt = dt, krylovdim = krylovdim, maxiter = maxiter, tol = tol)
        end 
    else
        @info "data exists! skip cal"
    end 

    if plot
        occplot(Par, filestr)
        curplot(Par, filestr)
    end 


    return nothing

end 