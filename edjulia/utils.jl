commutator(A, B) = A * B - B * A

anticommutator(A, B) = A * B + B * A

vectomat( vec ) = mapreduce( permutedims, vcat, vec)




function gen_file(top ; kwargs...)

    file = top * join([ string(key) * string(val) for (key, val) in kwargs], "_") * "/"

    return file

end 


function get_systag( sys)

    if typeof(sys) == TwoD
        L = sys.L
        dim = Int(sqrt(L))
        tag = "$(dim)x$(dim)"

    elseif typeof(sys) == SD
        X = sys.X
        dim = Int(X)
        tag = "SD$(dim)x$(dim)"
    end 

    return tag
end 


# function check_sparse(func, u0, p)

#     detector = SparseConnectivityTracer.TracerSparsityDetector()

#     @info "Real"

#     du0 = real.(u0)
#     jac_sparsity = ADTypes.jacobian_sparsity(
#     (du, u) -> (func(du, u, p, 0.0)), du0, u0, detector)

#     @info "Imag"
#     du0 = imag.(u0)
#     jac_sparsity = ADTypes.jacobian_sparsity(
#     (du, u) -> func(du, u, p, 0.0), du0, u0, detector)

# end 


