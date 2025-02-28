# function solve(h, Geo :: Geometry)

#     @assert h == transpose(h)
#     h = Hermitian(h)
#     eigen_h = @time eigen(h, 1:min(1, length(h)))

#     w = eigen_h.values
#     U = eigen_h.vectors

#     @show w
#     @show length(w)
#     @show w[1], w[end]

#     open( "ref/energy" * get_name(Geo), "w") do io
#         writedlm(io, w)
#     end

# end 



function solve(h, Geo :: Geometry, Par :: Particle)

    @assert h == transpose(h)
    h = Hermitian(h)
    #eigen_h = @time eigen(h, 1:min(1, length(h)))

    eigen_h, history = partialschur(h, nev =min(2, length(h)), tol=1e-7, which=:SR, restarts=1000)

    #display(eigen_h)

    @show w = real(eigen_h.eigenvalues)
    U = eigen_h.Q
    @show size(U)
    @show history
    #w, U = partialeigen(eigen_h)
    #w = eigen_h.values
    #U = eigen_h.vectors

    #@show w
    #@show length(w)
    #@show w[1], w[end]

    open( "ref/energy" * get_name(Geo) * get_name(Par), "w") do io
        writedlm(io, w)
    end

    return w, U

end 


