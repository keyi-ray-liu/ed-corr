function hopping!(basis_dict :: Dict, M :: Union{Matrix, SparseMatrixCSC}, ::Fermion, Geo :: Geometry; local_shift = 0)

    hopdict = Geo.hopdict

    for (basis :: Tuple, ind1 :: Int) in basis_dict
        for (from :: Int, tos :: Vector) in hopdict
            for to in tos
                

                from += local_shift
                to += local_shift

                @assert from < to

                # c^+ c only acts on larger
                if basis[from] < basis[to]
                    newbasis =  collect(basis)
                    newbasis[to], newbasis[from] = basis[from], basis[to]

                    ind2 = basis_dict[Tuple(newbasis)]
                    hop =  Geo.t * jordanwigner(Fermion(), basis, from, to)
                    M[ind1, ind2] += hop
                    M[ind2, ind1] += hop
                end 
            end 
        end 
    end 

    return M

end 


function hopping!(basis_dict :: Dict, M :: Union{Matrix, SparseMatrixCSC}, ::Electron, Geo :: Geometry; local_shift = 0, spin_shift = Geo.L)

    hopdict = Geo.hopdict
    for (basis :: Tuple, ind1 :: Int) in basis_dict
        for (from :: Int, tos :: Vector) in hopdict
            for to in tos
                
                from += local_shift
                to += local_shift

                @assert from < to

                #up
                if basis[from] < basis[to]
                    newbasis =  collect(basis)
                    newbasis[to], newbasis[from] = basis[from], basis[to]

                    ind2 = basis_dict[Tuple(newbasis)]
                    hop =  Geo.t * jordanwigner(Up(), basis, from, to, spin_shift)
                    M[ind1, ind2] += hop
                    M[ind2, ind1] += hop
                end 

                #dn
                if basis[from + spin_shift] < basis[to + spin_shift]
                    newbasis =  collect(basis)
                    newbasis[to + spin_shift], newbasis[from + spin_shift] = basis[from + spin_shift], basis[to + spin_shift]

                    ind2 = basis_dict[Tuple(newbasis)]
                    hop =  Geo.t * jordanwigner(Dn(), basis, from, to, spin_shift)
                    M[ind1, ind2] += hop
                    M[ind2, ind1] += hop
                end 


            end 
        end 
    end 

    return M

end 


function hopping!(basis_dict :: Dict, M :: Union{Matrix, SparseMatrixCSC}, Par::Electron, sd :: SD)

    # S, D hopping indices are faithful, no need to shift
    M = hopping!(basis_dict, M, Par, sd.S ; local_shift = 0, spin_shift = sd.L)
    M = hopping!(basis_dict, M, Par, sd.D ; local_shift = 0, spin_shift = sd.L)
    M = hopping!(basis_dict, M, Par, sd.A ; local_shift = sd.S.L + sd.D.L, spin_shift = sd.L)


end 




jordanwigner(::Fermion, basis, from, to) = (-1) ^ ( 1 + length(findall( >(1), basis[from + 1:to - 1] )))


function jordanwigner( ::Up, basis, from , to, L)


    up = length(findall( >(1), basis[(from ):(to - 1)] ))
    dn = length(findall( >(1), basis[(from + L):(to -1 + L)] ))

    total = up + dn
    
    #return  (-1)^ ( 1 + length(findall( x -> x == 2 || x == 4, basis[from + 1:to - 1] )))
    return  (-1) ^ (total)
end 

function jordanwigner( ::Dn, basis, from , to, L)


    up =  length(findall( >(1), basis[(from + 1):(to )] ))
    dn =  length(findall( >(1), basis[(from + 1 + L):(to + L )] ))
    
    total = up + dn
    #return  (-1)^ ( 1 + length(findall( x -> x == 2 || x == 4, basis[from + 1:to - 1] )))
    return  (-1) ^ (total)

end 