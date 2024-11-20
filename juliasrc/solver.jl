

vectomat( vec ) = mapreduce( permutedims, vcat, vec)
# BiasSolver class in Julia
struct BiasSolver
end

# Solve function equivalent in Julia
function solve!(::BiasSolver, h0::Matrix{Float64}, h::Matrix{Float64}, timestep::Float64, fin::Float64, N::Int, leftinds :: UnitRange, rightinds :: UnitRange; factor::Float64=1.0)
    # Eigen decomposition of matrices
    eigen_h0 = eigen(h0)
    U0 = eigen_h0.vectors
    eigen_h = eigen(h)
    w = eigen_h.values
    U = eigen_h.vectors

    total = length(w)
    steps = trunc(Int, fin / timestep)

    #E = zeros(ComplexF64, steps, total, total)

    # Populate E matrix
    # for t in 1:steps
    #     for i in 1:total
    #         for j in 1:total
    #             E[t, i, j] = exp(-1im * (w[i] - w[j]) * timestep * t)
    #         end
    #     end
    # end 

    # E = zeros(ComplexF64, steps, total, total)

    # for t in 1:steps
    #     E[t, :, :] = reshape(kron( exp.(1im * w * timestep * t), exp.(-1im * w * timestep * t)), (length(w), length(w)))
    # end 

    #@show count( f-> abs(f) > 1e-6, altE - E)

    E = zeros(ComplexF64, total, total)
    for i in 1:total
        for j in 1:total
            E[i, j] = exp(-1im * (w[i] - w[j]) * timestep)
        end
    end
    cur_E = copy(E)

    U0 = U0[:, 1:N]
    Ud = transpose(conj(U))
    U0d = transpose(conj(U0))

    UU = Ud * U0 * U0d * U
    #UL = Ud[:, leftinds]
    UL = transpose(Ud)[leftinds, :]
    #UR = U[rightinds, :]
    UR = transpose(U)[:, rightinds]


    CCs = Array{ComplexF64, 3}(undef, steps + 1, length(leftinds), length(rightinds))
    CCs[1, :, :] .= UL * ( UU.* ones(total, total)) * UR

    for i in 1:steps

        CC = UL * (UU .* cur_E) * UR
        CCs[i + 1, :, :] = CC
        cur_E = cur_E .* E
    end

    mem = Sys.maxrss()/GB
    @show "Max. RSS = $mem GB"

    return CCs, 0:steps * timestep

end


# Current calculation function
function current(sd::SD, CCs::Union{Matrix, Array{ComplexF64, 3}}; offset=0)
    currents = zeros(Float64, sd.arr * 2, size(CCs, 1))
    contact_source = sd.contact_source - offset
    contact_drain = sd.contact_drain - offset
    arr_source = sd.arr_source .- offset
    arr_drain = sd.arr_drain .- offset

    for i in 2:size(CCs, 1)  # Start from 2 as the initial state does not contribute to the current
        for j in 1:length(arr_source)
            arr = arr_source[j]
            cc = -2 * imag(CCs[i, contact_source, arr])
            currents[j, i] = cc
        end

        for k in 1:length(arr_drain)
            arr = arr_drain[k]
            cc = -2 * imag(CCs[i, arr, contact_drain]) 
            currents[k + length(arr_source), i] = cc
        end
    end

    return currents
end
