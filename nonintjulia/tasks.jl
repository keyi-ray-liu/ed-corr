
# SDbias function equivalent in Julia
function SDbias(;L::Int=128, N::Int=128, fin::Float64=128.0, biasSD::Float64=0.25, s::Float64=0.1, d::Float64=0.1, arr::Int=1, arrhop :: Float64 =1.0,single::Bool=true, plungerinit::Plunger=Plunger(arr), plungerdyna::Plunger=Pluger(arr), prefix::String="", mode="gs", ω=1.0)

    # time = now()
    timestep = 0.25
    key = ""
    #name = "L$(L)_N$(N)_dim$(arr)_biasSD$(biasSD)_$(plungerinit.name)_$(plungerdyna.name)_S$(s)_D$(d)_single$(single)"
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
        ("single", single),
        ("spacing", ω)
    ]

    name = join([ k == v ? k : (k * string(v)) for (k, v) in sort(namepairs)], "_")
    path = prefix * name 
    @show path

    # System initialization

    if mode == "gs"
        system = SD(L, N, arr, biasSD, plungerdyna.biases, -biasSD, s, d; arrhop=arrhop, single=single, ω=ω)
        INIT = SD(L, N, arr, 0.0, plungerinit.biases, 0.0, s, d; arrhop=arrhop, single=single, ω=ω)
    elseif mode == "left"
        system = SD(L, N, arr, 0.0, plungerdyna.biases, 0.0, s, d; arrhop=arrhop, single=single, ω=ω)
        INIT = SD(L, N, arr, -100.0 * ω, plungerinit.biases, 100.0 * ω, s, d; arrhop=arrhop, single=single, ω=ω)
    end 


    h0 = set_hij(INIT, path, "init")
    h = set_hij(system, path, "full")

    solver = BiasSolver()

    contact_source = system.contact_source
    contact_drain = system.contact_drain


    leftinds = rightinds = contact_source : contact_drain
    #leftinds = rightinds = 1:2 * L + arr^2

    CCs, _ = solve!(solver, h0, h, timestep, fin, N, leftinds, rightinds)
    @show size(CCs)

    currents = current(system, CCs; offset= leftinds[1] - 1)

    # Calculate occupancy of states in the array region
    occ = vectomat([abs.(diag(CCs[i, :, :])[2:end-1]) for i in axes(CCs, 1)])

    @show axes(occ)


    # Save results to files
    save_result(system, currents, occ, CCs, timestep, fin, path)

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
#     single = true
#     @show length(combs)

#     save_parameters(prefix; Ls=Ls, biasSDs=biasSDs, ss=ss, ds=ds, arrs=arrs, biasAs=biasAs, arrhops=arrhops)

#     Threads.@threads for (biasA, arr, arrhop, d, s, biasSD, L) ∈ combs
#     #for (biasA, arr, d, s, biasSD, L) ∈ combs

#         plungerinit = Plunger(arr; G1=biasAinit, G2=biasAinit, mid=biasAinit, name="biasAinit$(biasAinit)")
#         plungerdyna = Plunger(arr; G1=biasA, G2=biasA, mid=biasA, name="biasA$(biasA)")

#         @time SDbias(;L=L, N=L, fin=float(L), biasSD=biasSD, s=s, d=d, arr=arr, arrhop=arrhop, single=single, plungerinit=plungerinit, plungerdyna=plungerdyna, prefix=prefix)
#     end
# end


function Scalingtest()
    Ls = [128]
    #ss = [trunc(10.0^-2, digits=5)]
    #ds = [trunc(10.0^i, digits=5) for i in -3:-1]

    #ss = [ 2.0^i for i in -8:8]
    ss = [ 1/2, 1/32, 1/512]
    ds = ss
    arrs = [3]
    biasAs = [0.0]
    biasAinit= 1000.0
    
    ωs = [1.0, 1/16]
    prefix = "LeftSpacingScan/"

    # Generate combinations of parameters
    combs = collect(Iterators.product(biasAs, arrs, ωs, ds, ss,  Ls))
    single = true
    @show length(combs)

    save_parameters(prefix; Ls=Ls, ss=ss, ds=ds, arrs=arrs, biasAs=biasAs, ωs=ωs)

    Threads.@threads for (biasA, arr, ω, d, s, L) ∈ combs
    #for (biasA, arr, d, s, biasSD, L) ∈ combs

        plungerinit = Plunger(arr; G1=biasAinit, G2=biasAinit, mid=biasAinit, name="biasAinit$(biasAinit)")
        plungerdyna = Plunger(arr; G1=biasA, G2=biasA, mid=biasA, name="biasA$(biasA)")

        @time SDbias(;L=L, N=div(L, 2), fin=float(L),  s=s, d=d, arr=arr,  single=single, plungerinit=plungerinit, plungerdyna=plungerdyna, prefix=prefix, ω=ω, mode="left")
    end
end


function single()
    Ls = [128]
    biasSDs = [0.25]
    ss = [1/32]
    ds = [1/32]
    arrs = [ 1]
    biasAs = [0.0]
    #biasAinits = [1.0, 10.0, 100.0, 1000.0]
    biasAinits = [1000.0]

    # Generate combinations of parameters
    combs = collect(Iterators.product(biasAs, biasAinits, arrs, ds, ss, biasSDs, Ls))
    single = true

    for (biasA, biasAinit, arr, _, s, biasSD, L) in combs


        plungerinit = Plunger(arr; G1=biasAinit, G2=biasAinit, mid=biasAinit, name="biasAinit$(biasAinit)")
        plungerdyna = Plunger(arr; G1=biasA, G2=biasA, mid=biasA, name="biasA$(biasA)")

        @time SDbias(;L=L, N=64, fin=200.0, biasSD=biasSD, s=s, d=s, arr=arr, single=single, plungerinit=plungerinit, plungerdyna=plungerdyna, prefix="LeftSpacingScan/", mode="left")
    end
end



# function CSScan()
#     Ls = [64]
#     biasSDs = [1.0]
#     #ss = [trunc(10.0^-2, digits=5)]
#     #ds = [trunc(10.0^i, digits=5) for i in -3:-1]
#     ss = [0.1]
#     ds =[0.1]
#     arrs = [3]
#     G1s = -0.4:0.02:0.4
#     G2s = -0.4:0.02:0.4

#     biasAinit= 1000.0

#     prefix = "CSScan/"

#     # Generate combinations of parameters
#     combs = collect(Iterators.product(G1s, G2s, arrs, ds, ss, biasSDs, Ls))
#     single = true
#     @show length(combs)

#     save_parameters(prefix; Ls=Ls, biasSDs=biasSDs, ss=ss, ds=ds, arrs=arrs, G1s=G1s, G2s=G2s)

#     Threads.@threads for (G1, G2, arr, d, s, biasSD, L) ∈ combs
#     #for (biasA, arr, d, s, biasSD, L) ∈ combs

#         plungerinit = Plunger(arr; G1=biasAinit, G2=biasAinit, mid=biasAinit, name="biasAinit$(biasAinit)")
#         plungerdyna = Plunger(arr; G1=G1, G2=G2, mid=0.0, name="G1$(G1)_G2$(G2)")

#         @time SDbias(;L=L, N=L, fin=0.9*L, biasSD=biasSD, s=s, d=d, arr=arr, single=single, plungerinit=plungerinit, plungerdyna=plungerdyna, prefix=prefix)
#     end
# end
