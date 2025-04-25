
function gen_basis(::Conserved, Par::Fermion, Geo::Geometry)

    combs = combinations(1:Geo.L, Par.N)
    result = []
    
    for comb in combs
        arr = zeros(Int, Geo.L)
        arr[comb] .= 1
        push!(result, Tuple(arr))
    end 

    @show length(result)

    return result
end 


function gen_basis(::Conserved, Par::Electron, Geo::Geometry)

    L = Geo.L

    ups = combinations(1:L, Par.Nup)
    dns = combinations(L + 1: 2 *L, Par.Ndn)
    result = []
    
    
    for up in ups
        for dn in dns
            arr = zeros(Int, L * 2)
            arr[up] .= 1
            arr[dn] .= 1
            push!(result, Tuple(arr))
        end 
    end 

    @show length(result)

    return result   
end 


function gen_basis(::Not_conserved, ::Fermion, Geo:: Geometry)
    return Base.product([0:1 for _ in 1:Geo.L]...)
end 

function gen_basis(::Not_conserved, ::Electron, Geo:: Geometry)
    return Base.product([0:1 for _ in 1:(Geo.L * 2)]...)
end 


