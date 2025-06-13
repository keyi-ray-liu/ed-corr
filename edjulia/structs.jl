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
    t :: Float64
    function Line(L :: Int; t = -1.0, hopdict = nothing)
        #only hops to later in the chain

        if !(typeof(hopdict)  <: Dict)
            hopdict = Dict(
                i => [i + 1] for i in 1:L - 1
            )
        end 
        new(L, hopdict, t)
    end 
end 

struct TwoD <: Geometry
    X :: Int
    Y :: Int
    L :: Int
    hopdict :: Dict
    t :: Float64
    function TwoD(X, Y; t = -1.0)
        L = X * Y
        hopdict = Dict{Int, Array}(
            (x - 1) * Y + y => [] for x in 1:X for y in 1:Y
        )

        for x in 1:X
            for y in 1:Y - 1
                append!(hopdict[ (x -1) * Y  + y], [ (x - 1) * Y + y + 1])
            end 
        end 

        for x in 1:X - 1
            for y in 1:Y
                append!(hopdict[ (x - 1) * Y + y], [ x * Y + y])
            end 
        end 
        new(X, Y, L, hopdict, t)
    end 
end 




# the structure is S, D, A, hardcoded to have 1 site
struct SD <: Geometry
    S :: Line
    D :: Line
    A :: TwoD
    L :: Int
    function SD(X, Y; AS = [1], AD = [X * Y], scoup = 1.0, dcoup = 1.0)
        L = Ls + Ld + X * Y

        S = Line(Ls; hopdict = Dict( Ls => [ val + Ls + Ld for val in AS]), t = scoup)
        D = Line(Ld; hopdict = Dict( Ls + Ld => [ val + Ls + Ld for val in AD]), t = dcoup)
        A = TwoD(X, Y)
        new(S, D, A, L)
    end 
end 

get_name(Geo :: Line) = "Line" * string(Geo.L)
get_name(Geo:: TwoD) = "$(string(Geo.X))x$(string(Geo.Y))-single"
get_name(Geo::SD) = "SD$(get_name(Geo.A))"

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