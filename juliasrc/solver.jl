function solve(h, Geo :: Geometry)

    eigen_h = eigen(h)
    w = eigen_h.values
    U = eigen_h.vectors

    #@show w
    @show length(w)
    @show w[1], w[end]

    open( "ref/energy" * get_name(Geo), "w") do io
        writedlm(io, w)
    end

end 


