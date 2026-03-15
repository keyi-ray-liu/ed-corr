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

    elseif typeof(sys) == FRUSTRATION_4
        tag = "FRUSTRATION4"

    elseif typeof(sys) == FIVE
        tag = "FIVE"

    else
        error("$(sys): unknown type")
    end 

    return tag
end 


function load_corr(filestr)

    h5open( "$(filestr)CC.h5", "r") do f 

        fRU = f["REup"]
        fRD = f["REdn"]
        fIU = f["IMup"]
        fID = f["IMdn"]

        RU = read(fRU)
        RD = read(fRD)
        IU = read(fIU)
        ID = read(fID)

        return RU, RD, IU, ID

    end 

    

end 

function write_corr(filestr, CCups, CCdns; tag = "")

    T, Nx, Ny = size(CCups)


    h5open("$(filestr)CC$(tag).h5", "w") do f

        dRu = create_dataset(f, "REup", Float64, dataspace(T, Nx, Ny))
        dRd = create_dataset(f, "REdn", Float64, dataspace(T, Nx, Ny))
        dIu = create_dataset(f, "IMup", Float64, dataspace(T, Nx, Ny))
        dId = create_dataset(f, "IMdn", Float64, dataspace(T, Nx, Ny))

        write(dRu, real.(CCups))
        write(dRd, real.(CCdns))
        write(dIu, imag.(CCups))
        write(dId, imag.(CCdns))

    end 

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


