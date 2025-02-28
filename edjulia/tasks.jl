function task()

    Us = [2.0, 50, 100.0]
    Nup = 3

    # Geo = TwoD(3, 3)
    # Par = Fermion(4)
    # Coul = Coulomb(1.0, -1.0, 0.5, 0.5, 0.2)
    # task(Conserved(), Par, Geo, Coul)

    # Geo = TwoD(1, 1)
    # @show Geo
    # task(Not_conserved(), Geo)
    # Geo = Line(12)
    # Par = Fermion(6)
    # Coul = Coulomb(2.0, -1.0, 0.5, 0.5, 0.2)
    # task(Conserved(), Par, Geo, Coul)
    

    #taks(Conserved(), Geo, N)
    #task(Conserved(), L, N)

    res = []

    for U in Us

        Geo = TwoD(3, 3)
        Par = Electron(Nup, 4, U)
        Coul = Coulomb(0.0, -0.0, 0.5, 0.5, 0.0)
        prev = task(Conserved(), Par, Geo, Coul)[1]

        Geo = TwoD(3, 3)
        Par = Electron(Nup + 1, 4, U)
        Coul = Coulomb(0.0, -0.0, 0.5, 0.5, 0.0)
        plusone = task(Conserved(), Par, Geo, Coul)[1]

        append!(res, [U, plusone - prev])

    end 


    open("ref/$(Nup)plusone", "w") do io
        writedlm(io, res)
    end 
    


    # Geo = TwoD(3, 3)
    # Par = Electron(7, 2, 4.0)
    # Coul = Coulomb(0.0, -0.0, 0.5, 0.5, 0.0)
    # task(Conserved(), Par, Geo, Coul)


end 



task(qn::QN, Par :: Particle,  Geo::Geometry, Coul::Coulomb) = @time _task(qn::QN, Par :: Particle,  Geo::Geometry, Coul::Coulomb)



function _task(qn::QN, Par :: Particle,  Geo::Geometry, Coul::Coulomb)

    #Coul = Coulomb(2.0, -1.0, 0.5, 0.2)
    #Coul = Coulomb(0.0, -0.0, 0.5, 0.0)

    ham = gen_ham(qn, Par, Geo, Coul)
    w, U = solve(ham, Geo, Par)

    expectation(qn, Par, Geo, U, Occupation())

    return w

end 


