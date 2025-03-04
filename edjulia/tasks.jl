function GS()

    Us = 0:100
    Ns = [
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

    res = []
    energies = []

    for U in Us
        for (Nup, Ndn) in Ns
            Geo = TwoD(3, 3)
            Par = Electron(Nup, Ndn, U)
            Coul = Coulomb(0.0, -0.0, 0.5, 0.5, 0.0)
            (one, two )= GS(Conserved(), Par, Geo, Coul)

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



GS(qn::QN, Par :: Particle,  Geo::Geometry, Coul::Coulomb) = @time _GS(qn::QN, Par :: Particle,  Geo::Geometry, Coul::Coulomb)



function _GS(qn::QN, Par :: Particle,  Geo::Geometry, Coul::Coulomb)

    #Coul = Coulomb(2.0, -1.0, 0.5, 0.2)
    #Coul = Coulomb(0.0, -0.0, 0.5, 0.0)

    ham = gen_ham(qn, Par, Geo, Coul)
    w, U = _solve(ham, Geo, Par)

    expectation(qn, Par, Geo, U, Occupation())

    return w

end 




function time_evolve()

    Nup = 1 + 2
    Ndn = 1 + 2
    U = 4.0

    Geo = SD(1, 1, 3, 3)
    Par = Electron(Nup, Ndn, U)
    Coul = Coulomb(0.0, -0.0, 0.5, 0.5, 0.0)

    time_evolve(Conserved(), Par, Geo, Coul, LoadLeft(-1000, 0))


end 





function time_evolve(qn::QN, Par :: Particle,  Geo::SD, Coul::Coulomb, mode::LoadLeft)

    bias_init = Bias(vcat([ mode.init for _ in 1:Geo.S.L], [ 0 for _ in 1:Geo.L - Geo.S.L]))
    bias_dyna = Bias(vcat([0 for _ in 1:Geo.S.L + Geo.D.L], [mode.dyna for _ in 1:Geo.A.L]))

    timecontrol = TimeControl(100, 0.1)

    _time_evolve(qn, Geo, Par, Coul, bias_dyna, bias_init, timecontrol)

end 


