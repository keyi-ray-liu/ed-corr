function hubbard(basis_dict :: Dict, M :: Union{Matrix, SparseMatrixCSC}, par::Electron, Geo :: Geometry; shift = 0)

    L = Geo.L
    for (basis :: Tuple, ind :: Int) in basis_dict

        for s in 1:L 

            if basis[s] == basis[ s + L] == 2
                M[ind, ind] += par.U 
            end 
        end 
    end 

    return M

end 


function hubbard(basis_dict :: Dict, M :: Union{Matrix, SparseMatrixCSC}, par::Electron, sd :: SD)

    return hubbard(basis_dict, M, par, sd.A; shift = sd.S.L)


end 