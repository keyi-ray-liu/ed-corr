function coulomb!(basis_dict :: Dict, M :: Union{Matrix, SparseMatrixCSC}, par::Particle, Coul :: Coulomb, Geo :: Geometry)

    for (basis :: Tuple, ind :: Int) in basis_dict

        M[ind, ind] += ee(basis, par, Coul, Geo)
    end 

    return M

end 


function ee(basis, ::Fermion, Coul::Coulomb, Geo :: Geometry)

    val = 0
    nonzero = findall( >(1), basis)
    comb = combinations(nonzero, 2)


    for (ind1, ind2) in comb

        dist = distance(ind1, ind2, Geo)
        r = dist + Coul.ζ_ee
        factor = dist == 1 ? 1 - Coul.exch : 1
        val += Coul.ee * factor / r
    end 


    for ind in eachindex(basis)
        for ind2 in nonzero

            if ind != ind2

                r = distance(ind, ind2, Geo) + Coul.ζ_ne
                val += Coul.ne / r
            end 
        end 
    end     

    return val

end 


function ee(basis, ::Electron, Coul::Coulomb, Geo::Geometry)

    val = 0
    nonzero = findall( x -> x > 1, basis)
    comb = combinations(nonzero, 2)


    for (ind1, ind2) in comb

        multiplier = div(basis[ind1], 2) * div(basis[ind2], 2)

        dist = distance(ind1, ind2, Geo)
        r = dist + Coul.ζ_ee
        factor = dist == 1 ? 1 - Coul.exch : 1
        val += Coul.ee * factor * multiplier / r
    end 


    for ind in eachindex(basis)
        for ind2 in nonzero

            if ind != ind2

                multiplier = 2 * div(basis[ind2], 2)
                r = distance(ind, ind2, Geo) + Coul.ζ_ne
                val += Coul.ne * multiplier / r
            end 
        end 
    end     

    return val

end 



distance(ind1, ind2, ::Line) = abs(ind1 - ind2)
function distance(ind1, ind2, Geo::TwoD)

    x1 = (ind1 - 1) ÷ Geo.Y
    y1 = (ind1 - 1) % Geo.Y

    x2 = (ind2 - 1) ÷ Geo.Y
    y2 = (ind2 - 1) % Geo.Y

    return sqrt( (x1 - x2)^2 + (y1 - y2)^2)


end 