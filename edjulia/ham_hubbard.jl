# local shift refers to the # of sites we shift due to preceding geo (reservoir for example), spin shift is the total geo size

function hubbard(basis_dict :: Dict, M :: Union{Matrix, SparseMatrixCSC}, par::Electron, Geo :: Geometry; local_shift = 0, spin_shift = Geo.L)

    L = Geo.L
    for (basis :: Tuple, ind :: Int) in basis_dict

        for s in 1 + local_shift : L + local_shift

            if basis[s] == basis[ s + spin_shift ] == 2
                M[ind, ind] += par.U 
            end 
        end 
    end 

    return M

end 


function hubbard(basis_dict :: Dict, M :: Union{Matrix, SparseMatrixCSC}, par::Electron, sd :: SD)

    return hubbard(basis_dict, M, par, sd.A; local_shift = sd.S.L + sd.D.L, spin_shift = sd.L)


end 