
let
    
    #Gs = -1000.0:20.0:1000.0

    Gs = -10.0:0.25:10.0

    open("commandlines", "w") do io
        for (G1, G2) in collect(Base.product(Gs, Gs))
            write(io, "julia main.jl GG $(G1) $(G2) \n")
        end 
    end 

end 