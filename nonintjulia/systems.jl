const GB = 2^30




struct Connector
    arr_source :: Array
    source_scaling :: Array
    arr_drain :: Array
    drain_scaling :: Array
    name :: String
    function Connector(; arr_source = [4], arr_drain = [6], source_scaling = [1.0], drain_scaling = [1.0], name = "46")
        new(arr_source, source_scaling, arr_drain, drain_scaling, name )
    end
end 

struct Plunger
    name :: String
    biases :: Vector{Float64}
    function Plunger(arr :: Int ; G1=0.0, G2=0.0, mid = 0.0, override=nothing, name=nothing)

        if isnothing(name)
            error("Plunger needs a name")
        end 

        if isnothing(override)
            biases = ones( arr^2) .* mid
            for first in 1:arr
                biases[first] = G1
            end 

            for last in arr^2 - arr + 1: arr^2
                biases[last] = G2
            end 
        else
            if typeof(override) != Vector{Float64}
                error("wrong type override")
            end 

            if len(override) != arr ^ 2 
                error("wrong dim of override")
            end 

            biases = override
        end 
        new(name, biases)
    end 
end 


# SD class in Julia
struct SD
    L::Int
    N::Int
    arr::Int
    left_bias::Float64
    arr_bias::Vector{Float64}
    right_bias::Float64
    source_coupling::Float64
    drain_coupling::Float64
    t0::Float64
    arrhop ::Float64
    contact_source::Int
    contact_drain::Int
    arr_source::Vector{Int}
    source_scaling::Vector{Float64}
    arr_drain::Vector{Int}
    drain_scaling::Vector{Float64}

    function SD(L::Int, N::Int, arr::Int, left_bias::Float64, arr_bias::Vector{Float64}, right_bias::Float64, source_coupling::Float64, drain_coupling::Float64; con :: Connector = Connector(), arrhop ::Float64 = 1.0, t0::Float64=-1.0)
        contact_source = L 
        contact_drain = L + arr^2 + 1
        arr_source, source_scaling, arr_drain, drain_scaling = setup_connections(con, contact_source, contact_drain)
        new(L, N, arr, left_bias, arr_bias, right_bias, source_coupling, drain_coupling, t0, arrhop, contact_source, contact_drain, arr_source, source_scaling, arr_drain, drain_scaling)
    end

    # Function to set up source and drain connections based on single flag and array size
    function setup_connections(  con :: Connector, contact_source::Int, contact_drain::Int)


        arr_source = con.arr_source .+ contact_source   
        arr_drain = con.arr_drain .+ contact_source
        return arr_source, con.source_scaling, arr_drain, con.drain_scaling
    end
end

# Function to set up the Hamiltonian matrix
function set_hij(sd::SD, path::String, mod::String; savehij::Bool=false)
    t0 = sd.t0
    arrhop = sd.arrhop
    L = sd.L
    arr = sd.arr
    arr_size = arr^2
    left = sd.left_bias
    cent = sd.arr_bias
    right = sd.right_bias
    source_coupling = sd.source_coupling
    drain_coupling = sd.drain_coupling
    contact_source = sd.contact_source
    contact_drain = sd.contact_drain
    arr_source = sd.arr_source
    arr_drain = sd.arr_drain
    source_scaling = sd.source_scaling
    drain_scaling = sd.drain_scaling

    total = 2 * L + arr_size
    M = zeros(Float64, total, total)

    # Left section
    for i in 1:(L - 1)
        M[i + 1, i] = M[i, i + 1] = t0
        M[i, i] = left
    end
    M[L, L] = left

    # Coupling between sections
    for i in 1:length(arr_source)
        site = arr_source[i]
        M[contact_source, site] = M[site, contact_source] = source_coupling * source_scaling[i]
    end
    for i in 1:length(arr_drain)
        site = arr_drain[i]
        M[contact_drain, site] = M[site, contact_drain] = drain_coupling *  drain_scaling[i]
    end

    # Right section
    for i in (L + arr_size + 1):(total - 1)
        M[i + 1, i] = M[i, i + 1] = t0
        M[i, i] = right
    end
    M[total, total] = right

    # Central array section
    for i in 1:arr
        for j in 1:arr
            idx = (i - 1) * arr + j + L
            M[idx, idx] += cent[idx - L]
            if i < arr
                M[idx, idx + arr] = arrhop
                M[idx + arr, idx] = arrhop
            end
            if j < arr
                M[idx, idx + 1] = arrhop
                M[idx + 1, idx] = arrhop
            end
        end
    end

    if M != transpose(conj(M))
        throw(DomainError("Matrix M is not Hermitian!"))
    end

    if savehij
        open(path * "/M" * mod * ".txt", "w") do f
            writedlm(f, M)
        end
    end

    return M
end

# Function to initialize states
function set_init(sd::SD)
    init = zeros(ComplexF64, sd.L * 2 + sd.arr^2)
    init[1:sd.N] .= 1.0 # Set the first N elements
    return init
end



# Helper function to save results
function save_result(::SD, currents::Array{Float64, 2}, occs::Array{Float64, 2}, CCs::Array{ComplexF64, 3}, time, path::String)
    mkpath(path)

    open(path * "/times", "w") do f
        writedlm(f, time)
    end
    open(path * "/currentSD", "w") do f
        writedlm(f, round.(currents; sigdigits=5))
    end
    open(path * "/occ", "w") do f
        writedlm(f, round.(occs; sigdigits=5))
    end
end