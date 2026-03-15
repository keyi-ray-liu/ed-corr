abstract type QN end 
abstract type Particle end
abstract type Geometry end
abstract type Spin end
abstract type Observable end
abstract type Parameter end
abstract type drivingmode end
abstract type DiagMethod end
abstract type OpenSystem end 

struct Up <: Spin end
struct Dn <: Spin end


struct full_diag <: DiagMethod end 
struct arnoldi <: DiagMethod end 

struct LoadLeft <: drivingmode 
    init :: Number
    dyna :: Number
end

struct InjDep <: OpenSystem
    source_site :: Int
    drain_site :: Int
    γ_inj_source :: Float64
    γ_dep_source :: Float64
    γ_inj_drain :: Float64
    γ_dep_drain :: Float64
end 


struct LoadBoth <: drivingmode

    init :: Number
    dyna :: Number
end 


struct Conserved <: QN
end 

struct Not_conserved <: QN
end 

struct Fermion <: Particle 
    N :: Int
    function Fermion(N=0)
        new(N)
    end 
end 

struct Electron <: Particle 
    Nup :: Int
    Ndn :: Int
    U :: Float64
    function Electron(;Nup = 0, Ndn = 0, U = 0)
        new(Nup, Ndn, U)
    end 
    Electron(Nup, Ndn, U) = Electron(; Nup = Nup, Ndn = Ndn, U = U)
end

get_name(Par :: Fermion) = "N$(Par.N)"
get_name(Par:: Electron) = "$(Par.Nup)Up$(Par.Ndn)DnU$(Par.U)"

struct Line <: Geometry 
    L :: Int
    hopdict :: Dict
    function Line(L :: Int; t = -1.0, hopdict = nothing)
        #only hops to later in the chain

        if !(typeof(hopdict)  <: Dict)
            hopdict = Dict(
                i => [(t, i + 1)] for i in 1:L - 1
            )
        end 
        new(L, hopdict)
    end 
end 

struct TwoD <: Geometry
    X :: Int
    Y :: Int
    L :: Int
    hopdict :: Dict
    function TwoD(X, Y; t = -1.0)
        L = X * Y
        hopdict = Dict{Int, Array}(
            (x - 1) * Y + y => [] for x in 1:X for y in 1:Y
        )

        for x in 1:X
            for y in 1:Y - 1
                append!(hopdict[ (x -1) * Y  + y], [ (t, (x - 1) * Y + y + 1)])
            end 
        end 

        for x in 1:X - 1
            for y in 1:Y
                append!(hopdict[ (x - 1) * Y + y], [ (t, x * Y + y)])
            end 
        end 
        new(X, Y, L, hopdict)
    end 
end 

"""1, is contact, 4 is drain, 5 is center, \n
   2
 1-5-4
   3
"""
struct FIVE <: Geometry
    L :: Int
    hopdict :: Dict
    function FIVE(; t = -1.0, tc = -1.0)
        L = 5
        hopdict = Dict{Int, Array}(
            1 => [ (t, 2), (t, 3), (tc, 5)],
            2 => [ (t, 4), (tc, 5)],
            3 => [ (t, 4), (tc, 5)],
            4 => [(tc, 5)],
            5 => []

        )
        new( L, hopdict)
    end 
end 


"""1, is contact, 4 is drain, but 1 has connection to 4"""
struct FRUSTRATION_4 <: Geometry
    L :: Int
    hopdict :: Dict
    function FRUSTRATION_4(; t = -1.0)
        L = 5
        hopdict = Dict{Int, Array}(
            1 => [ (t, 2), (t, 3), (t, 4)],
            2 => [ (t, 4)],
            3 => [ (t, 4)],
            4 => []

        )
        new(L, hopdict)
    end 
end 




# the structure is S, D, A
struct SD <: Geometry
    L :: Int
    X :: Int
    Y :: Int
    RES :: Int
    hopdict :: Dict
    scoup :: Float64
    dcoup :: Float64
    function SD(X, Y; scoup = -1.0, dcoup = -1.0, RES = 1, ω = -1.0)
        L = 2 * RES  + X * Y
        A = TwoD(X, Y)

        hopdict = A.hopdict
        
        newhop = Dict()

        for (key, value) in hopdict
            newhop[key + 2 * RES] = [ (t, v + 2 * RES) for (t, v) in value]
        end 

        # connect

        # last site of source
        # first site of drain
        newhop[RES ] = [ (scoup, 2 * RES + 1)]
        newhop[RES + 1] = [ (dcoup, L)]

        # internal

        for nxt in 1:(RES - 1)
            #source
            newhop[nxt] = [ (ω, nxt + 1)]

            #drain
            newhop[nxt + RES] = [ (ω, nxt + 1 + res)]
        end 
        
        @show newhop
        new(L, X, Y, RES, newhop, scoup, dcoup)
    end 
end 

get_name(Geo :: Line) = "Line" * string(Geo.L)
get_name(Geo:: TwoD) = "$(string(Geo.X))x$(string(Geo.Y))-single"
get_name(Geo::SD) = "SD$(string(Geo.X))x$(string(Geo.Y))res$(string(Geo.L))"

struct Coulomb  <: Parameter
    ee :: Float64
    ne :: Float64
    ζ_ee :: Float64
    ζ_ne :: Float64
    exch :: Float64
end 


struct Bias <: Parameter
    val :: Array{Number}
end 

struct TimeControl <: Parameter
    fin :: Number
    dt :: Number
end 





RefCoul =  Coulomb(2.0, -1.0, 0.5, 0.5, 0.2)
StrongCoul = Coulomb(200.0, -100.0, 0.5, 0.5, 0.2)
ZeroCoul = Coulomb(0.0, 0.0, 0.5, 0.5, 0.2)