function hopping!(basis_dict :: Dict, M :: Matrix, ::Fermion, Geo :: Geometry)

    hopdict = Geo.hopdict

    for (basis :: Tuple, ind1 :: Int) in basis_dict
        for (from :: Int, tos :: Vector) in hopdict
            for to in tos

                if basis[from] != basis[to]
                    newbasis =  collect(basis)
                    newbasis[to], newbasis[from] = basis[from], basis[to]

                    ind2 = basis_dict[Tuple(newbasis)]
                    M[ind1, ind2] += -1
                end 
            end 
        end 
    end 

    return M

end 