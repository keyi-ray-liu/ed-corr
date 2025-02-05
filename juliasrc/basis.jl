
function gen_basis(::Conserved, Sp::Fermion, Geo::Geometry)

    combs = combinations(1:Geo.L, Sp.N)
    result = []
    
    for comb in combs
        arr = zeros(Int, Geo.L)
        arr[comb] .= 1
        push!(result, Tuple(arr))
    end 

    return result
end 


function gen_basis(::Not_conserved, ::Fermion, Geo:: Geometry)
    return Base.product([0:1 for _ in 1:Geo.L]...)
end 

function gen_basis(::Not_conserved, ::Electron, Geo:: Geometry)
    return Base.product([0:3 for _ in 1:Geo.L]...)
end 


