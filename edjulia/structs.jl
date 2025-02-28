abstract type QN end 
abstract type Particle end
abstract type Geometry end
abstract type Spin end
abstract type Observable end


struct Up <: Spin end
struct Dn <: Spin end





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
end

get_name(Par :: Fermion) = "N$(Par.N)"
get_name(Par:: Electron) = "$(Par.Nup)Up$(Par.Ndn)DnU$(Par.U)"

struct Line <: Geometry 
    L :: Int
    hopdict :: Dict
    t :: Float64
    function Line(L :: Int; t = -1.0)
        #only hops to later in the chain
        hopdict = Dict(
            i => [i + 1] for i in 1:L - 1
        )
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

get_name(Geo :: Line) = "Line" * string(Geo.L)
get_name(Geo:: TwoD) = "$(string(Geo.X))x$(string(Geo.Y))-single"



struct SD <: Geometry
    S :: Line
    D :: Line
    A :: TwoD
    L :: Int
    AS :: Array
    AD :: Array
    function SD(S, D, A; AS = [4], AD = [6])
        L = S.L + D.L + A.L
        new(S, D, A, L, AS, AD)
    end 
end 

struct Coulomb 
    ee :: Float64
    ne :: Float64
    ζ_ee :: Float64
    ζ_ne :: Float64
    exch :: Float64
end 
