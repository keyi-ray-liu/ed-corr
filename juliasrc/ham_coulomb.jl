function coulomb!(basis_dict :: Dict, M :: Matrix, sp::Fermion, Coul :: Coulomb)

    for (basis :: Tuple, ind :: Int) in basis_dict

        M[ind, ind] += ee(basis, sp, Coul)
    end 

    return M

end 


function ee(inds, ::Fermion, Coul::Coulomb)

    val = 0
    nonzero = findall( x -> x > 0, inds)
    comb = combinations(nonzero, 2)


    for (ind1, ind2) in comb

        dist = abs(ind1 - ind2)
        r = dist + Coul.ζ
        factor = dist == 1 ? 1 - Coul.exch : 1
        val += Coul.ee * factor / r
    end 


    for ind in eachindex(inds)
        for ind2 in nonzero

            if ind != ind2

                r = abs(ind - ind2) + Coul.ζ
                val += Coul.ne / r
            end 
        end 
    end     

    return val

end 

