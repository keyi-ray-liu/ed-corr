function _hopping(basis_dict :: Dict, M :: SparseMatrixCSC, ::Fermion, hopdict)


    for (basis :: Tuple, ind1 :: Int) in basis_dict
        for (from :: Int, tos :: Vector) in hopdict
            for (t, to) in tos
                

                @assert from < to

                # c^+ c only acts on larger
                if basis[from] < basis[to]
                    newbasis =  collect(basis)
                    newbasis[to], newbasis[from] = basis[from], basis[to]

                    ind2 = basis_dict[Tuple(newbasis)]
                    hop =  t * jw(Fermion(), basis, from, to)
                    M[ind1, ind2] += hop
                    M[ind2, ind1] += hop
                end 
            end 
        end 
    end 

    return M

end 


function hopping!(basis_dict :: Dict, M :: SparseMatrixCSC, ::Fermion, Geo :: Geometry)

    hopdict = Geo.hopdict

    M = _hopping(basis_dict, M, Fermion(), hopdict)

    return M

end 

# generates the hopping matrix
function _hopping(basis_dict :: Dict, M :: SparseMatrixCSC, ::Electron, hopdict::Dict,   spin_shift::Number ; up=true, dn = true, herm = true )

    for (basis :: Tuple, ind1 :: Int) in basis_dict
        for (from :: Int, tos :: Vector) in hopdict
            for (t, to) in tos

                if up
                    if basis[from] < basis[to]
                        newbasis =  collect(basis)
                        newbasis[to], newbasis[from] = basis[from], basis[to]

                        ind2 = basis_dict[Tuple(newbasis)]
                        hop =  t* jw(Up(), basis, from, to, spin_shift)
                        M[ind1, ind2] += hop

                        if herm
                            M[ind2, ind1] += hop
                        end 
                    end 
                end 

                #dn
                if dn
                    if basis[from + spin_shift] < basis[to + spin_shift]
                        newbasis =  collect(basis)
                        newbasis[to + spin_shift], newbasis[from + spin_shift] = basis[from + spin_shift], basis[to + spin_shift]

                        ind2 = basis_dict[Tuple(newbasis)]
                        hop =  t * jw(Dn(), basis, from, to, spin_shift)
                        M[ind1, ind2] += hop

                        if herm
                            M[ind2, ind1] += hop
                        end 
                    end
                end 


            end 
        end 
    end 

    return M

end 


function hopping!(basis_dict :: Dict, M :: SparseMatrixCSC, Par::Electron, Geo :: Geometry; spin_shift = Geo.L)

    hopdict = Geo.hopdict

    M = _hopping(basis_dict, M, Par, hopdict,  spin_shift)

    return M

end 





# ---------- new --------------



jw(::Fermion, basis, site) = (-1) ^ ( 1 + length(findall( >(0), basis[1:site - 1] )))

function jw( ::Up, basis, site, L)


    up = length(findall( >(0), basis[ 1 : site - 1] ))
    dn = length(findall( >(0), basis[(1 + L):( site - 1 + L)] ))

    total = up + dn
    
    #return  (-1)^ ( 1 + length(findall( x -> x == 2 || x == 4, basis[from + 1:to - 1] )))
    return  (-1) ^ (total)
end 

function jw( ::Dn, basis, site, L)


    up = length(findall( >(0), basis[ 1 : site ] ))
    dn = length(findall( >(0), basis[(1 + L):( site - 1 + L)] ))

    total = up + dn
    
    #return  (-1)^ ( 1 + length(findall( x -> x == 2 || x == 4, basis[from + 1:to - 1] )))
    return  (-1) ^ (total)
end 

jw( ::Fermion, basis, from, to) = jw(Fermion(), basis, from) * jw(Fermion(), basis, to)
jw( Sp, basis, from, to, L) = jw(Sp, basis, from, L) * jw(Sp, basis, to, L)



# hop_amp(td :: TwoD, :: Int, :: Int) = td.t
# hop_amp(ln :: Line, :: Int, :: Int) = ln.t


# function hopping_direct!(basis_dict, M,  :: Fermion, Geo)

#     hopdict = Geo.hopdict

#     for (from, tos) in hopdict
#         for to in tos
#             M += hop_amp(Geo, from, to) * ( cdag(basis_dict, Geo, from) * c(basis_dict, Geo, to) + cdag(basis_dict, Geo, to) * c(basis_dict, Geo, from) ) 
#         end 
#     end 

#     return M

# end 




# function hopping_direct!(basis_dict, M,  :: Electron, Geo)

#     hopdict = Geo.hopdict

#     for (from, tos) in hopdict
#         for to in tos
#             M += hop_amp(Geo, from, to)* ( cdagup(basis_dict, Geo, from) * cup(basis_dict, Geo, to) + cdagup(basis_dict, Geo, to) * cup(basis_dict, Geo, from) )
#             M += hop_amp(Geo, from, to) * ( cdagdn(basis_dict, Geo, from) * cdn(basis_dict, Geo, to) + cdagdn(basis_dict, Geo, to) * cdn(basis_dict, Geo, from) )
#         end 
#     end 

#     return M

# end 