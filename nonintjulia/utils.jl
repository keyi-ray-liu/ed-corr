function save_parameters(prefix; args...)

    for (key, val) in args
        
        name = String(key)
        open(prefix * name, "w") do io
            writedlm(io, collect(val))
        end 

    end 
end 


function slide(M)
    fig = Figure()

    ax = Axis(fig[1, 1])

    sl_x = Slider(fig[2, 1], range = 1:128, startvalue = 1)

    val = lift(sl_x.value) do x
        M[:, x]
    end

    scatter!(val, color = :red, markersize = 20)
    #limits!(ax, 0, 10, 0, 10)

    display(fig)
end 


function gradient(x, y)
    itp = linear_interpolation(x, y)

    val = Interpolations.gradient.(Ref(itp), x)

    val = [ v[1] for v in val]
    return val
end 