function task()



    Geo = TwoD(2, 2)
    @show Geo
    task(Not_conserved(), Geo)

    Geo = TwoD(3, 3)
    @show Geo
    task(Not_conserved(), Geo)

    Geo = TwoD(1, 1)
    @show Geo
    task(Not_conserved(), Geo)
    #Geo = Line(L)
    #taks(Conserved(), Geo, N)
    #task(Conserved(), L, N)

end 




function task(::Not_conserved, Geo::Geometry)

    N = 0
    type = Fermion(N) 
    Coul = Coulomb(2.0, -1.0, 0.5, 0.2)


    ham = gen_ham(Not_conserved(), type, Geo, Coul)
    solve(ham, Geo)

    return nothing

end 


function task(::Conserved, Geo::Geometry, N)

    type = Fermion(N) 
    Coul = Coulomb(2.0, -1.0, 0.5, 0.2)


    ham = gen_ham(Conserved(), type, Geo, Coul)
    solve(ham, Geo)

    return nothing

end 