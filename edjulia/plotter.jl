function occplot(::Electron, filestr)
    
    plt = plot() 

    Nups = open( "$(filestr)occup", "r") do io
        readdlm(io)
    end

    Ndns = open( "$(filestr)occdn", "r") do io
        readdlm(io )
    end
    
    t= open( "$(filestr)time", "r") do io
        readdlm(io)
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
    #display(plt)


end 

function occplot(::Fermion)

    plt = plot()
    
    for (i, obs) in enumerate(eachcol(Ns))
        plot!(plt, sol.t, obs; label = "n$(i)", xlabel = "time", title = "", ylabel = "⟨n̂⟩")
    end 

    savefig(plt, "test/testocc.pdf")
    #display(plt)


end 



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

    curs = open( "$(filestr)current", "r") do io
        readdlm(io)
    end

    
    t = open( "$(filestr)time", "r") do io
        readdlm(io)
    end


    for (i, obs) in enumerate(eachcol(curs))
        plot!(plt, t, obs; label = "cur $(i)", xlabel = "time", title = "ED ODE test Current", ylabel = "I")
    end 

    savefig(plt, "$(filestr)testcur.pdf")
    display(plt)

end 