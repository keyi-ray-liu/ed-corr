
function gen_basis(::Conserved, Par::Fermion, Geo::Geometry)

    combs = combinations(1:Geo.L, Par.N)
    result = []
    
    for comb in combs
        arr = ones(Int, Geo.L)
        arr[comb] .= 2
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
            arr = ones(Int, L * 2)
            arr[up] .= 2
            arr[dn] .= 2
            push!(result, Tuple(arr))
        end 
    end 

    @show length(result)

    return result   
end 


function gen_basis(::Not_conserved, ::Fermion, Geo:: Geometry)
    return Base.product([1:2 for _ in 1:Geo.L]...)
end 

function gen_basis(::Not_conserved, ::Electron, Geo:: Geometry)
    return Base.product([1:2 for _ in 1:(Geo.L * 2)]...)
end 


