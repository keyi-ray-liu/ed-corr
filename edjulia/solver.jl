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



function _solve(h, Geo :: Geometry, Par :: Particle; nev=min(2, length(h)))

    @assert h == transpose(h)
    h = Hermitian(h)
    #eigen_h = @time eigen(h, 1:min(1, length(h)))

    eigen_h, history = partialschur(h, nev =nev, tol=1e-7, which=:SR, restarts=1000)

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


function _time_evolve(qn::QN, Geo::Geometry, Par::Particle, Coul::Coulomb, bias_dyna::Bias, bias_init::Bias, timecontrol :: TimeControl; obs=[])

    @show Geo
    h_dyna = gen_ham(qn, Par, Geo, Coul, bias_dyna )
    h_init = gen_ham(qn, Par, Geo, Coul, bias_init)

    _, v = _solve(h_init, Geo, Par; nev=1)
    w, U = _solve(h_dyna, Geo, Par; nev=40)

    # we find the coefficients of the all the overlaps, then explicitly time_evolve

    coeffs = vec(conj(transpose(v)) * U)

    @show size(coeffs)
    @show sum( abs2, coeffs)


    times = []
    occups = []
    occdns = []

    for t in 0:timecontrol.dt : timecontrol.fin
        append!(times, t)
        newcoeff = coeffs .* exp.( im * w * t)

        V = U * newcoeff
        V = reshape(V, (size(V)..., 1))

        @show size(newcoeff), size(V), typeof(V)
        occup, occdn = expectation(qn, Par, Geo, V, Occupation())

        append!(occups, occup)
        append!(occdns, occdn)
    
    end 

    filestr = "ref/$(get_name(Geo))$(get_name(Par))/"

    try
        mkdir( filestr )
    catch
    end

    open( "$(filestr)TEoccup", "w") do io
        writedlm(io, occups)
    end

    open( "$(filestr)TEoccdn", "w") do io
        writedlm(io, occdns)
    end

    open( "$(filestr)time", "w") do io
        writedlm(io, times)
    end


end 

