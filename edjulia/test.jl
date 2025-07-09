



function time_evolve()

    # Nup = 4
    # Ndn = 5
    # U = 100.0
    # gap = -1.79161798

    Nup = 2
    Ndn = 2
    U = 4.0
    gap = 0.0

    # Nup = 1
    # Ndn = 1
    # U = 4.0
    # gap = 0.0

    Geo_init = SD(1, 1, 3, 3; scoup = 0, dcoup = 0)
    Geo_dyna = SD(1, 1, 3, 3; scoup = -1/8, dcoup = -1/8)
    Par = Electron(Nup, Ndn, U)
    Coul = Coulomb(0.0, -0.0, 0.5, 0.5, 0.0)

    filestr = "ref/$(get_name(Geo_dyna))$(get_name(Par))bias$(gap)/"

    time_evolve(Conserved(), Par, Geo_init, Geo_dyna, Coul, LoadLeft(-1000, gap), filestr)
    #time_evolve(Conserved(), Par, Geo_init, Geo_dyna, Coul, LoadBoth(1000, 0))

end 





function time_evolve(qn::QN, Par :: Particle,  Geo_init::SD, Geo_dyna::SD, Coul::Coulomb, mode::LoadLeft, filestr::String)

    bias_init = Bias(vcat([ mode.init for _ in 1:Geo_init.S.L], [ -mode.init for _ in 1:Geo_init.D.L], [ 0 for _ in 1:Geo_init.A.L]))
    bias_dyna = Bias(vcat([0 for _ in 1:Geo_dyna.S.L + Geo_dyna.D.L], [mode.dyna for _ in 1:Geo_dyna.A.L]))

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
    Geo = TwoD(3,3)
    Par = Fermion(4)
    Coul = Coulomb(0.0, 0.0, 0.5, 0.5, 0.2)
    bias = Bias([0 for _ in 1:Geo.L])
    qn = Conserved()

    M1 = gen_ham(qn, Par, Geo, Coul, bias)
    M2 = gen_ham_direct(qn, Par, Geo, Coul, bias)


    isapprox(M1, M2)

end 





function markov_test()


    #Geo = TwoD(2, 2)
    
    Geo = SD(2, 2; scoup = -0.1, dcoup = -0.1 )
    systag = get_systag(Geo)
    #Par = Fermion() 

    Par = Electron(; 
    U = 4.0
    )

    Coul = Coulomb(
    2.0,
    -1.0, 
    0.5, 0.5, 0.2)#[0.1, 1.0, 10.0, 100.0]

    γs = [0.5, 0.01]

    G1 = 0.0
    G2 = 0.0

    devicebias = 0.0
    bias = Bias( [0, 0, 0, G1, G2, 0] .+ devicebias)

    states = [
        (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
        (1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0),
    ]

    for γ in γs

        injdeps = [
            InjDep(1, 2, γ, 0, 0, γ)
        ]

        for injdep in injdeps
            for state in states
                top = "/Users/knl20/Desktop/PROJECT_SD/Markovian_ED/$(systag)/"
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


                @show filestr

                if !ispath(filestr * "time")
                    ρ = gen_ρ(Not_conserved(), Par, Geo; state = state)
                    @time odesolve(Not_conserved(), Par, Geo, Coul, bias, ρ , injdep; filestr = filestr, start = 0, fin = 5, chunks = 1)

                    occplot(Par, filestr)
                    curplot(Par, filestr)
                end 
            end 
        end 
    end 

    return nothing
end 




