function GS()

    Us = 0:100
    Ns = [
        (0, 1),
        (1, 1),
        (1, 2),
        (0, 3),
        (1, 6),
        (2, 7),
        (3, 4),
        (1, 7),
        (2, 6),
        (3, 5),
        (4, 4)
    ]

    # Geo = TwoD(3, 3)
    # Par = Fermion(4)
    # Coul = Coulomb(1.0, -1.0, 0.5, 0.5, 0.2)
    # GS(Conserved(), Par, Geo, Coul)

    # Geo = TwoD(1, 1)
    # @show Geo
    # GS(Not_conserved(), Geo)
    # Geo = Line(12)
    # Par = Fermion(6)
    # Coul = Coulomb(2.0, -1.0, 0.5, 0.5, 0.2)
    # GS(Conserved(), Par, Geo, Coul)
    

    #taks(Conserved(), Geo, N)
    #GS(Conserved(), L, N)

    #res = []
    energies = []

    for U in Us
        for (Nup, Ndn) in Ns
            Geo = TwoD(3, 3)
            Par = Electron(Nup, Ndn, U)
            Coul = Coulomb(0.0, -0.0, 0.5, 0.5, 0.0)
            bias = Bias( [0 for _ in 1:Geo.L])
            (one, two )= GS(Conserved(), Par, Geo, Coul, bias)

            #append!(res, [U, plusone - prev])
            append!(energies, [[Nup, Ndn, U, one, two]])

        end 
    end 


    # open("ref/$(Nup)plusone", "w") do io
    #     writedlm(io, res)
    # end 
    
    open("ref/energyscan", "w") do io
        writedlm(io, energies)
    end 
    


    # Geo = TwoD(3, 3)
    # Par = Electron(7, 2, 4.0)
    # Coul = Coulomb(0.0, -0.0, 0.5, 0.5, 0.0)
    # GS(Conserved(), Par, Geo, Coul)


end 



GS(qn::QN, Par :: Particle,  Geo::Geometry, Coul::Coulomb, bias::Bias) = @time _GS(qn, Par,Geo, Coul, bias)



function _GS(qn::QN, Par :: Particle,  Geo::Geometry, Coul::Coulomb, bias :: Bias)

    #Coul = Coulomb(2.0, -1.0, 0.5, 0.2)
    #Coul = Coulomb(0.0, -0.0, 0.5, 0.0)

    ham = gen_ham(qn, Par, Geo, Coul, bias)
    w, U = _solve(ham, Geo, Par)

    expectation(qn, Par, Geo, U, Occupation())

    return w

end 




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

    timecontrol = TimeControl(1000, 0.1)

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







function scan()


    re = readdlm("ref/energyscan")
    refdict = Dict()


    for (Nup, Ndn, U, gs, _) in eachrow(re)

        refdict[ (Int(Nup), Int(Ndn), U) ] = gs 

        if Ndn != Nup
            refdict[ (Int(Ndn), Int(Nup), U)] = gs
        end 
    end 



    Us = 0:10:100

    for U in Us
        refdict[(0, 0, U)] = 0.0
    end 

    Ns = [
        (0, 0),
        (1, 1),
        (2, 0),
        #(1, 6),
        #(3, 4)
    ]

    multipliers = -2:0.2:2
    scs = [-1.0, -0.5, -0.1, -0.01, -0.001]
    dcs = [-1.0, -0.5, -0.1, -0.01, -0.001]

    for (Nup, Ndn) in Ns
        for U in Us
            Threads.@threads for multiplier in multipliers
                for sc in scs
                    for dc in dcs
                        Geo_init = SD(1, 1, 3, 3; scoup = 0, dcoup = 0)
                        Geo_dyna = SD(1, 1, 3, 3; scoup = sc, dcoup = dc)
                        Par = Electron(Nup + 1, Ndn + 1, U)
                        Coul = Coulomb(0.0, -0.0, 0.5, 0.5, 0.0)

                        gap = (refdict[ (Nup + 1, Ndn, U)] - refdict[ (Nup, Ndn, U)]) * multiplier

                        filestr = "scan/$(get_name(Geo_dyna))_$(get_name(Par))_1emultiplier$(multiplier)_sc$(sc)_dc$(dc)/"
                    
                        time_evolve(Conserved(), Par, Geo_init, Geo_dyna, Coul, LoadLeft(-1000, gap), filestr)
                    end 
                end 
            end 
        end 
    end 
    

end 