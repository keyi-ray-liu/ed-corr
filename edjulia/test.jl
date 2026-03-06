



function time_evolve_test()

    # Nup = 4
    # Ndn = 5
    # U = 100.0
    # gap = -1.79161798

    Nup = 1
    Ndn = 1
    U = 4.0
    gap = 0.0

    # Nup = 1
    # Ndn = 1
    # U = 4.0
    # gap = 0.0

    Geo_init = SD(2, 2; scoup = 0, dcoup = 0, res = 20)
    Geo_dyna = SD(2, 2; scoup = -1/8, dcoup = -1/8, res = 20)
    Par = Electron(Nup, Ndn, U)
    Coul = Coulomb(0.0, -0.0, 0.5, 0.5, 0.0)

    filestr = "ref/$(get_name(Geo_dyna))$(get_name(Par))bias$(gap)/"

    time_evolve(Conserved(), Par, Geo_init, Geo_dyna, Coul, LoadLeft(-1000, gap), filestr)
    #time_evolve(Conserved(), Par, Geo_init, Geo_dyna, Coul, LoadBoth(1000, 0))

end 





function time_evolve(qn::QN, Par :: Particle,  Geo_init::SD, Geo_dyna::SD, Coul::Coulomb, mode::LoadLeft, filestr::String)


    res = div(Geo_init.L - Geo_init.X * Geo_init.Y, 2)
    arr = Geo_init.X * Geo_init.Y

    bias_init = Bias(vcat([ mode.init for _ in 1:res], [ -mode.init for _ in 1:res], [ 0 for _ in 1:arr]))
    bias_dyna = Bias(vcat([0 for _ in 1:2 * res], [mode.dyna for _ in 1:arr]))

    timecontrol = TimeControl(1000, 0.25)

    @show bias_init
    @show bias_dyna

    _time_evolve(qn, Geo_init, Geo_dyna, Par, Coul, bias_dyna, bias_init, timecontrol, [Occupation(), Current()], filestr)

end 


# function time_evolve(qn::QN, Par :: Particle,  Geo_init::SD, Geo_dyna::SD, Coul::Coulomb, mode::LoadBoth)

#     bias_init = Bias(vcat([ mode.init for _ in 1:Geo_init.S.L + Geo_init.D.L], [ 0 for _ in 1:Geo_init.A.L]))
#     bias_dyna = Bias(vcat([0 for _ in 1:Geo_dyna.S.L + Geo_dyna.D.L], [mode.dyna for _ in 1:Geo_dyna.A.L]))

#     @show bias_init
#     timecontrol = TimeControl(100, 0.1)

#     _time_evolve(qn, Geo_init, Geo_dyna, Par, Coul, bias_dyna, bias_init, timecontrol, [Occupation()])

# end 

function setup_test()

    sd = SD(2, 2; scoup = -0.1, dcoup = -0.1)

    @show sd.hopdict

end 






function test()

    Coul = Coulomb(2.0, -1.0, 0.5, 0.5, 0.2)

    Geo = Line(12)
    Par = Fermion(6)
    

    #Coul = Coulomb(0.0, 0.0, 0.1, 0.1, 0.0)
    bias = Bias( [0 for _ in 1:Geo.L])
    val = GS(Conserved(), Par, Geo, Coul, bias; nev=1)

    @show val

    Geo = TwoD(3, 3)
    Par = Fermion(4)
    
    #Coul = Coulomb(0.0, 0.0, 0.1, 0.1, 0.0)
    bias = Bias( [0 for _ in 1:Geo.L])
    val = GS(Conserved(), Par, Geo, Coul, bias; nev=1)

    @show val


    Geo = TwoD(3, 3)
    Par = Electron(3, 4, 4.0)
    
    #Coul = Coulomb(0.0, 0.0, 0.1, 0.1, 0.0)
    bias = Bias( [0 for _ in 1:Geo.L])
    val = GS(Conserved(), Par, Geo, Coul, bias; nev=1)

    @show val



    Geo = TwoD(3, 3)
    Coul = Coulomb(0.0, -0.0, 0.5, 0.5, 0.2)
    Par = Electron(3, 4, 100.0)
    
    #Coul = Coulomb(0.0, 0.0, 0.1, 0.1, 0.0)
    bias = Bias( [0 for _ in 1:Geo.L])
    val = GS(Conserved(), Par, Geo, Coul, bias; nev=1)

    @show val
    # open("ref/$(Nup)plusone", "w") do io
    #     writedlm(io, res)
    # end 
    

end 




function check_ham_construction()


    Geo = TwoD(2, 4)
    Par = Electron(; U =4.0)
    Coul = Coulomb(2.0, 2.0, 0.5, 0.5, 0.2)
    bias = Bias([0 for _ in 1:Geo.L])
    qn = Not_conserved()

    M1 = gen_ham(qn, Par, Geo, Coul, bias)
    # M2 = gen_ham_direct(qn, Par, Geo, Coul, bias)


    # isapprox(M1, M2)

    open("/Users/knl20/Desktop/Code/ED/ham/matrix$(Geo.X)x$(Geo.Y)", "w") do IO
        writedlm(IO, M1)
    end 

end 




function markov_test_SD()


    #Geo = TwoD(2, 2)
    Geo = SD(2, 2; scoup = -0.02, dcoup = -0.02)
    systag = get_systag(Geo)
    #Par = Fermion() 

    Par = Electron(; U = 1.0)

    Coul = ZeroCoul

    G1 = G2 = 0
    devicebias = 0
    γ = 0.5
    bias = Bias([0, 0, 0, 0, 0, 0])

    state = (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,)
    injdep = InjDep(1, 2, γ, 0.0 , 0.0, γ)

    top = "/home/keyi-liu/Desktop/Code/Markovian/Mar2test/$(systag)/"
    filestr = gen_file(top; 
        U = Par.U,
        Coul = Coul.ee,
        injs = injdep.γ_inj_source,
        deps = injdep.γ_dep_source,
        injd = injdep.γ_inj_drain,
        depd = injdep.γ_dep_drain,
        GOne = G1,
        GTwo = G2,
        devicebias = devicebias,
        state = join(state, "")
    )

    @show state
    @show filestr

    if !ispath(filestr * "time")
        ρ = gen_ρ(Not_conserved(), Par, Geo; state = state)
        @time odesolve(Not_conserved(), Par, Geo, Coul, bias, ρ , injdep; filestr = filestr, start = 0, fin = 10, chunks = 1)

    else
        @info "data exists! skip cal"
    end 


    occplot(Par, filestr)
    #curplot(Par, filestr)

    return nothing
end 




function markov_test()


    Geo = TwoD(2, 2)
    systag = get_systag(Geo)
    #Par = Fermion() 

    Par = Electron(; U = 1.0)

    Coul = ZeroCoul

    G1 = G2 = 0
    devicebias = 0
    γ = 0.5
    bias = Bias([0, 0, 0, 0])

    state = (0, 0, 0, 0, 0, 0, 0, 0,)
    injdep = InjDep(1, 4, γ, 0.0 , 0.0, γ)

    top = "test/$(systag)/" #"/home/keyi-liu/Desktop/Code/Markovian/Mar2test/$(systag)/"
    filestr = gen_file(top; 
        U = Par.U,
        Coul = Coul.ee,
        injs = injdep.γ_inj_source,
        deps = injdep.γ_dep_source,
        injd = injdep.γ_inj_drain,
        depd = injdep.γ_dep_drain,
        GOne = G1,
        GTwo = G2,
        devicebias = devicebias,
        state = join(state, "")
    )

    @show state
    @show filestr

    if !ispath(filestr * "time")
        ρ = gen_ρ(Not_conserved(), Par, Geo; state = state)
        @time odesolve(Not_conserved(), Par, Geo, Coul, bias, ρ , injdep; filestr = filestr, start = 0, fin = 10, chunks = 1, backendstr = "generic")

    else
        @info "data exists! skip cal"
    end 


    occplot(Par, filestr)
    #curplot(Par, filestr)

    return nothing
end 





# function speedtest()

#     dim = 200

#     A = rand((dim, dim))
#     L = 
# end 

