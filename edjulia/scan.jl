function _GG(G1, G2)
    Geo = TwoD(2, 2)
    Par = Electron(;U = 4.0)
    Coul = Coulomb(2.0, -1.0, 0.5, 0.5, 0.2)
    γ = 2.0

    G1 = trunc(G1; sigdigits = 5)
    G2 = trunc(G2; sigdigits = 5)

    filestr = "GGscan/g$(γ)_GOne$(G1)_GTwo$(G2)_U$(Par.U)_Coul$(Coul.ee)/"

    if !ispath(filestr * "time")
    
        #Par = Electron(;U = 4.0)
        bias = Bias( [0, G1, G2, 0])
        ρ = gen_ρ(Not_conserved(), Par, Geo)
        @time odesolve(Not_conserved(), Par, Geo, Coul, bias, ρ ,InjDep(1, 4, γ, 0.0, 0.0, γ); filestr = filestr, start = 0, fin = 500, chunks = 1)
    

    end 
end 


function GG()
    
    
    #Gs = -1000.0:20.0:1000.0

    Gs = -10.0:1.0:10.0
    Threads.@threads for (G1, G2) in collect(Base.product(Gs, Gs))
        _GG(G1, G2)
    end 

end 


function GG(G1, G2)
    _GG(G1, G2)
end 




function gamma_scan(; Geo = nothing)


    
    γs = 10.0 .^ (-2:0.1:3.0) #[0.1, 1.0, 10.0, 100.0]
    Us = [1.0, 10.0, 100.0]
    tfin = 500
    dt = 500

    krylovdim = 90
    maxiter = 500
    tol = 1e-11

    for U in Us

        U = trunc(U; sigdigits = 5)

        if isnothing(Geo)
            Geo = TwoD(2, 2)
        end 
        #Geo = SD(2, 2; scoup = -0.02, dcoup = -0.02)
        systag = get_systag(Geo)
        solver = "exp"
        #Par = Fermion() 
        Threads.@threads for γ in γs

            γ = trunc(γ; sigdigits = 5)
            Par = Electron(; U = U)

            Coul = ZeroCoul

            G1 = G2 = 0
            devicebias = 0
            bias = Bias([0 for _ in 1:Geo.L])

            state = Tuple([0 for _ in 1:Geo.L * 2])
            injdep = InjDep(1, 4, γ, 0.0 , 0.0, γ)

            top = "/home/keyi-liu/Desktop/Code/Markovian/Mar13EXP/"
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


            #occplot(Par, filestr)
            #curplot(Par, filestr)
        end 
    end 

    return nothing
end 