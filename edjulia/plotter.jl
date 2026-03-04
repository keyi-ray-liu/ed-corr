function occplot(::Electron, filestr)
    
    plt = plot() 

    # Nups = open( "$(filestr)occup", "r") do io
    #     readdlm(io)
    # end

    # Ndns = open( "$(filestr)occdn", "r") do io
    #     readdlm(io )
    # end
    
    t= open( "$(filestr)time", "r") do io
        readdlm(io)
    end

    RU, RD, _, _ = load_corr(filestr)

    T, D, _ = size(RU)

    Nups = zeros(T, D)
    Ndns = zeros(T, D)

    for a in 1:T, i in 1:D
        Nups[a, i] = RU[a, i, i]
        Ndns[a, i] = RD[a, i, i]
    end 


    for (i, obs) in enumerate(eachcol(Nups))
        
        if i != 3
            plot!(plt, t, obs; label = "n↑$(i)", xlabel = "time", title = "", ylabel = "⟨n̂⟩")
        else
            scatter!(plt, t, obs; label = "n↑$(i)", xlabel = "time", title = "", ylabel = "⟨n̂⟩", marker= :circle)
        end 
    end 

    for (i, obs) in enumerate(eachcol(Ndns))
        plot!(plt, t, obs; label = "n↓$(i)", xlabel = "time", title = "ED ODE test ⟨n⟩", ylabel = "⟨n̂⟩", linestyle = :dash)
    end 

    savefig(plt, "$(filestr)testocc.pdf")
    display(plt)


end 

# function occplot(::Fermion)

#     plt = plot()
    
#     for (i, obs) in enumerate(eachcol(Ns))
#         plot!(plt, sol.t, obs; label = "n$(i)", xlabel = "time", title = "", ylabel = "⟨n̂⟩")
#     end 

#     savefig(plt, "test/testocc.pdf")
#     #display(plt)


# end 



function curplot(::Fermion)

    plt = plot()
    

    for (i, obs) in enumerate(eachcol(curs))
        plot!(plt, sol.t, obs; label = "cur $(i)", xlabel = "time", title = "", ylabel = "I")
    end 

    savefig(plt, "test/testcur.pdf")
    display(plt)

end 


function curplot(::Electron, filestr)
    
    plt = plot()

    _, _, IU, ID = load_corr(filestr)

    T, _, _ = size(IU)

    # sites 1, 2, 3, 4

    curs = zeros(T, 4)

    for a in 1:T
        curs[a, 1] = 2 * (IU[a, 1, 2] + IU[a, 1, 3] + IU[a, 1, 4] )
        curs[a, 2] = 2 * (IU[a, 1, 4] + IU[a, 2, 4] + IU[a, 3, 4] )
        curs[a, 3] = 2 * (ID[a, 1, 2] + ID[a, 1, 3] + ID[a, 1, 4] )
        curs[a, 4] = 2 * (ID[a, 1, 4] + ID[a, 2, 4] + ID[a, 3, 4] )

    end 
    
    t = open( "$(filestr)time", "r") do io
        readdlm(io)
    end

    for (i, obs) in enumerate(eachcol(curs))
        plot!(plt, t, obs; label = "cur $(i)", xlabel = "time", title = "ED ODE test Current \n $(filestr)", ylabel = "I")
    end 

    savefig(plt, "$(filestr)testcur.pdf")
    display(plt)

end 