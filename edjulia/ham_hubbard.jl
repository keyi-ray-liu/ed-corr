# local shift refers to the # of sites we shift due to preceding geo (reservoir for example), spin shift is the total geo size

function hubbard!(basis_dict :: Dict, M :: SparseMatrixCSC, par::Electron, Geo :: Geometry)

    L = Geo.L
    for site in 1:L
        M += hubbardU(Geo, par, site) * nupdn(basis_dict, Geo, site)
    end 

    return M
end 



hubbardU(::Geometry, Par :: Electron, site) = Par.U
hubbardU( sd :: SD, Par :: Electron, site) = site > 2 * sd.RES ? Par.U : 0