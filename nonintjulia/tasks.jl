
# SDbias function equivalent in Julia
function SDbias(;L::Int=128, N::Int=128, fin::Float64=128.0, biasSD::Float64=0.25, s::Float64=0.1, d::Float64=0.1, arr::Int=1, arrhop :: Float64 =-1.0, con::Connector = Connector(), plungerinit::Plunger=Plunger(arr), plungerdyna::Plunger=Pluger(arr), prefix::String="", mode="gs", t0=1.0, manual_current=false)

    
    # time = now()
    timestep = 0.25
    key = ""
    namepairs = [
        ("L", L),
        ("N", N),
        ("biasSD", biasSD),
        ("arr", arr),
        ("arrhop", arrhop),
        (plungerinit.name, plungerinit.name),
        (plungerdyna.name, plungerdyna.name),
        ("s", s),
        ("d", d),
        ("con", con.name),
        ("spacing", t0),
        ("mode", mode),
    ]

    name = join([ k == v ? k : (k * string(v)) for (k, v) in sort(namepairs)], "_")
    path = prefix * name 
    @show path

    # System initialization

    if mode == "gs"
        system = SD(L, N, arr, biasSD, plungerdyna.biases, -biasSD, s, d; arrhop=arrhop, con=con, t0=t0)
        INIT = SD(L, N, arr, 0.0, plungerinit.biases, 0.0, s, d; arrhop=arrhop, con=con, t0=t0)
    elseif mode == "left"
        system = SD(L, N, arr, 0.0, plungerdyna.biases, 0.0, s, d; arrhop=arrhop, con=con, t0=t0)
        INIT = SD(L, N, arr, - 1e5 * abs(t0), plungerinit.biases, 1000.0 * t0, s, d; arrhop=arrhop, con=con, t0=t0)
    end 


    h0 = set_hij(INIT, path, "init")
    h = set_hij(system, path, "full")




    solver = BiasSolver()

    contact_source = system.contact_source
    contact_drain = system.contact_drain

    if !manual_current
        leftinds = rightinds = contact_source : contact_drain
    else
        leftinds = rightinds = 1:(2 * L + arr^2)
    end 

    CCs, _ = solve!(solver, h0, h, timestep, fin, N, leftinds, rightinds)
    @show size(CCs)

    currents = current(system, CCs; offset= leftinds[1] - 1)
    occ = vectomat([abs.(diag(CCs[i, :, :])) for i in axes(CCs, 1)])

    time = collect(0:timestep:fin)
    # Save results to files
    save_result(system, currents, occ, CCs, time,  path)

    if manual_current

        leftocc = sum(occ[:, 1: contact_source], dims=2)
        rightocc = sum(occ[:, contact_drain:end], dims=2)
        cur = -gradient(collect(0:timestep:fin), vec(leftocc))
        cur[1] = 0.0

        cur2 = gradient(collect(0:timestep:fin), vec(rightocc))
        cur2[1] = 0.0


        plt = plot()
        plot!(plt, time, cur ; label ="manual")
        plot!(plt, time, currents[:, 1]; label = "auto")

        plot!(plt, time, cur2 ; label ="manualdrain")
        plot!(plt, time, currents[:, 2]; label = "autodrain")
        title!(plt, "s$(s)d$(d)t0$(t0)")
        savefig(plt, path * "/comp.png")
        
        open(path * "/manualcurrent", "w") do f
            writedlm(f, round.(cur; sigdigits=5))
        end

        open(path * "/leftocc", "w") do f
            writedlm(f, round.(leftocc; sigdigits=5))
        end

    end 


    # avg_time = (now() - time) / Millisecond(1000) / nworkers()
    # @info "avg time ", avg_time
    return nothing

end




# function SmallScan()
#     Ls = [256, 512]
#     #ss = [trunc(10.0^-2, digits=5)]
#     #ds = [trunc(10.0^i, digits=5) for i in -3:-1]

#     ss = [ 4.0^-i for i in 1:5]
#     ds = [4.0^-i for i in 1:5]
#     arrs = [1, 2, 3, 4]
#     #biasAs = 0.0:0.1:3.0
#     biasAs = [0.0]
#     biasAinit= 1000.0

#     prefix = "SDScan/"

#     # Generate combinations of parameters
#     combs = collect(Iterators.product(biasAs, arrs, arrhops, ds, ss, biasSDs, Ls))
#     @show length(combs)

#     save_parameters(prefix; Ls=Ls, biasSDs=biasSDs, ss=ss, ds=ds, arrs=arrs, biasAs=biasAs, arrhops=arrhops)

#     Threads.@threads for (biasA, arr, arrhop, d, s, biasSD, L) ∈ combs
#     #for (biasA, arr, d, s, biasSD, L) ∈ combs

#         plungerinit = Plunger(arr; G1=biasAinit, G2=biasAinit, mid=biasAinit, name="biasAinit$(biasAinit)")
#         plungerdyna = Plunger(arr; G1=biasA, G2=biasA, mid=biasA, name="biasA$(biasA)")

#         @time SDbias(;L=L, N=L, fin=float(L), biasSD=biasSD, s=s, d=d, arr=arr, arrhop=arrhop, plungerinit=plungerinit, plungerdyna=plungerdyna, prefix=prefix)
#     end
# end



function manual()
    Ls = [128]

    ss = [-0.5, -0.1]
    ds = ss
    arrs = [3]
    biasAs = [0.0]
    biasAinit= 1000.0
    t0s = [-0.1, -1.0, -10.0]
    
    prefix = "manualtest/"

    if !isdir(prefix)
        mkdir(prefix)
    end 

    # Generate combinations of parameters
    combs = collect(Iterators.product(biasAs, arrs, t0s, ds, ss,  Int.(Ls)))
    @show length(combs)

    con = Connector(; arr_source = [1], arr_drain = [9], name = "19")

    save_parameters(prefix; Ls=Ls, ss=ss, ds=ds, arrs=arrs, biasAs=biasAs, t0s=t0s)

    for (biasA, arr, t0, d, s, L) ∈ combs
    #for (biasA, arr, d, s, biasSD, L) ∈ combs

        plungerinit = Plunger(arr; G1=biasAinit, G2=biasAinit, mid=biasAinit, name="biasAinit$(biasAinit)")
        plungerdyna = Plunger(arr; G1=biasA, G2=biasA, mid=biasA, name="biasA$(biasA)")

        @time SDbias(;L=L, N=L, fin=600.0, biasSD=0.0, s=s, d=d, arr=arr, con=con, plungerinit=plungerinit, plungerdyna=plungerdyna, prefix=prefix, t0=t0, mode="left", manual_current = true)
        #@time SDbias(;L=L, N=L, fin=600.0, biasSD=0.5 * abs(t0) , s=s, d=d, arr=arr, con=con, plungerinit=plungerinit, plungerdyna=plungerdyna, prefix=prefix, t0=t0, mode="gs")
    end
end




function LengthScale()
    Ls = 2.0 .^ (0:9)

    ss = [-100.0, -0.5, -0.1, -0.01]
    ds = ss
    arrs = [3]
    biasAs = [0.0]
    biasAinit= 1000.0
    t0s = [-0.1, -1.0, -10.0]
    
    prefix = "LengthScale/"

    if !isdir(prefix)
        mkdir(prefix)
    end 

    # Generate combinations of parameters
    combs = collect(Iterators.product(biasAs, arrs, t0s, ds, ss,  Int.(Ls)))
    @show length(combs)

    con = Connector(; arr_source = [1], arr_drain = [9], name = "19")

    save_parameters(prefix; Ls=Ls, ss=ss, ds=ds, arrs=arrs, biasAs=biasAs, t0s=t0s)

    Threads.@threads for (biasA, arr, t0, d, s, L) ∈ combs
    #for (biasA, arr, d, s, biasSD, L) ∈ combs

        plungerinit = Plunger(arr; G1=biasAinit, G2=biasAinit, mid=biasAinit, name="biasAinit$(biasAinit)")
        plungerdyna = Plunger(arr; G1=biasA, G2=biasA, mid=biasA, name="biasA$(biasA)")

        @time SDbias(;L=L, N=L, fin=600.0, biasSD=0.0, s=s, d=d, arr=arr, con=con, plungerinit=plungerinit, plungerdyna=plungerdyna, prefix=prefix, t0=t0, mode="left")
        @time SDbias(;L=L, N=L, fin=600.0, biasSD=0.5 * abs(t0) , s=s, d=d, arr=arr, con=con, plungerinit=plungerinit, plungerdyna=plungerdyna, prefix=prefix, t0=t0, mode="gs")
    end
end