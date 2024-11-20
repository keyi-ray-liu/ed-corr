
# SDbias function equivalent in Julia
function SDbias(;L::Int=128, N::Int=128, fin::Float64=128.0, bias::Float64=0.25, s::Float64=0.1, d::Float64=0.1, arr::Int=1, single::Bool=true, mu_init::Float64=1000.0, mu_te::Float64=0.0, prefix::String="")

    # time = now()
    timestep = 0.125
    key = "SDbias"
    name = "L$(L)_N$(N)_$(arr)x$(arr)_bias$(bias)_muinit$(mu_init)_mute$(mu_te)_S$(s)_D$(d)_single$(single)"
    path = prefix * name * "_" * key

    # System initialization
    system = SD(L, N, arr, bias, mu_te, -bias, s, d; single=single)
    INIT = SD(L, N, arr, 0.0, mu_init, 0.0, s, d; single=single)

    h0 = set_hij(INIT, path, "init")
    h = set_hij(system, path, "full")

    solver = BiasSolver()

    contact_source = system.contact_source
    contact_drain = system.contact_drain

    leftinds = rightinds = contact_source : contact_drain

    CCs, _ = solve!(solver, h0, h, timestep, fin, N, leftinds, rightinds)

    currents = current(system, CCs; offset= contact_source - 1)

    # Calculate occupancy of states in the array region
    occ = vectomat([abs.(diag(CCs[i, :, :])[2:end-1]) for i in axes(CCs, 1)])


    # Save results to files
    save_result(system, currents, occ, CCs, timestep, fin, path)

    # avg_time = (now() - time) / Millisecond(1000) / nworkers()
    # @info "avg time ", avg_time
    return nothing

end


# SDscan function to run simulations with multiple parameter combinations
function SDscan()
    Ls = [2^i for i in 6:10]
    biases = [2.0^i for i in -3:3]
    ss = [trunc(10.0^i, digits=5) for i in -3:1]
    ds = [trunc(10.0^i, digits=5) for i in -3:1]
    arrs = [1, 2, 3, 4]
    mu_tes = [0; 2.0.^(-3:4)]

    @show mu_tes

    # Generate combinations of parameters
    combs = collect(Iterators.product(mu_tes, arrs, ds, ss, biases, Ls))
    single = true
    @show length(combs)

    open("parameters", "w") do io
        writedlm(io, combs)
    end 

    Threads.@threads for (mu_te, arr, d, s, bias, L) ∈ combs
    #for (mu_te, arr, d, s, bias, L) ∈ combs
        @time SDbias(;L=L, N=L, fin=0.9*L, bias=bias, s=s, d=d, arr=arr, single=single, mu_init=1000.0, mu_te=mu_te, prefix="SDscanjlnew/")
    end
end


function SmallScan()
    Ls = [64]
    biases = 0.0:0.1:1.0
    ss = [trunc(10.0^-2, digits=5)]
    ds = [trunc(10.0^i, digits=5) for i in -4:-2]
    arrs = [1, 2, 3, 4]
    mu_tes = 0.0:0.1:1.0

    @show mu_tes

    # Generate combinations of parameters
    combs = collect(Iterators.product(mu_tes, arrs, ds, ss, biases, Ls))
    single = true
    @show length(combs)

    open("parameters", "w") do io
        writedlm(io, combs)
    end 

    Threads.@threads for (mu_te, arr, d, s, bias, L) ∈ combs
    #for (mu_te, arr, d, s, bias, L) ∈ combs
        @time SDbias(;L=L, N=L, fin=0.9*L, bias=bias, s=s, d=d, arr=arr, single=single, mu_init=1000.0, mu_te=mu_te, prefix="SDscanjlsmall/")
    end
end



function single()
    Ls = [128]
    biases = [0.25]
    ss = [0.1]
    ds = [0.1]
    arrs = [ 1, 2, 3]
    mu_tes = [1.0]

    # Generate combinations of parameters
    combs = collect(Iterators.product(mu_tes, arrs, ds, ss, biases, Ls))
    single = true

    for (mu_te, arr, d, s, bias, L) in combs
        @time SDbias(;L=L, N=L, fin=200.0, bias=bias, s=s, d=d, arr=arr, single=single, mu_init=1.0, mu_te=mu_te, prefix="")
    end
end