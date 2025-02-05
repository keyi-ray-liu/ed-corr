abstract type QN end 
abstract type Spin end
abstract type Geometry end

struct Conserved <: QN
end 

struct Not_conserved <: QN
end 

struct Fermion <: Spin 
    N :: Int
end 

struct Electron <: Spin 
    Nup :: Int
    Ndn :: Int
    U :: Float64
end

struct Line <: Geometry 
    L :: Int
    hopdict :: Dict
    function Line(L :: Int)
        #only hops to later in the chain
        hopdict = Dict(
            i => [i + 1] for i in 1:L - 1
        )
        new(L, hopdict)
    end 
end 

struct TwoD <: Geometry
    X :: Int
    Y :: Int
    L :: Int
    hopdict :: Dict
    function TwoD(X, Y)
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
        new(X, Y, L, hopdict)
    end 
end 

get_name(Geo :: Line) = "Line" * string(Geo.L)
get_name(Geo:: TwoD) = "$(string(Geo.X))x$(string(Geo.Y))-single"


struct Coulomb 
    ee :: Float64
    ne :: Float64
    ζ :: Float64
    exch :: Float64
end 