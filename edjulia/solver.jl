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

    @time eigen_h, history = partialschur(h, nev =nev, tol=1e-7, which=:SR, restarts=1000)

    @show w = real(eigen_h.eigenvalues)
    U = eigen_h.Q
    @show size(U)
    @show history

    # w, U, info = eigsolve(h,  nev, :SR)

    # @show info
    # @show size(w), size(U)


    filestr = "ref/$(get_name(Geo))$(get_name(Par))/"

    try
        mkdir( filestr )
    catch
    end

    open( "$(filestr)energy" , "w") do io
        writedlm(io, w)
    end

    return w, U

end 


function _time_evolve(qn::QN, Geo_init::Geometry, Geo_dyna::Geometry, Par::Particle, Coul::Coulomb, bias_dyna::Bias, bias_init::Bias, timecontrol :: TimeControl, obs :: Array{T} where T <: Observable, filestr::String)

    @show Geo_init
    @show Geo_dyna


    times = collect(0:timecontrol.dt : timecontrol.fin)

    h_dyna = gen_ham(qn, Par, Geo_dyna, Coul, bias_dyna )
    h_init = gen_ham(qn, Par, Geo_init, Coul, bias_init)

    _, v = _solve(h_init, Geo_init, Par; nev=1)
    w, U = _solve(h_dyna, Geo_dyna, Par; nev=min(800, size(v, 1) - 2))

    # we find the coefficients of the all the overlaps, then explicitly time_evolve

    coeffs = vec(conj(transpose(v)) * U)
    E = exp.(im * w * times')

    @show size(coeffs)
    @show size(E)
    @show sum( abs2, coeffs)


    newcoeffs = coeffs' .* E'
    @show size(newcoeffs)

    V = U * newcoeffs'
    @show size(V)


    for ob in obs
        expectation(qn, Par, Geo_dyna, V, ob; filestr =filestr)
    end 


    try
        mkdir( filestr )
    catch
    end

    open( "$(filestr)times", "w") do io
        writedlm(io, times)
    end


end 

