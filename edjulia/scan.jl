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



