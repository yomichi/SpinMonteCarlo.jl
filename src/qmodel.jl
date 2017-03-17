@enum(LocalOperatorType,
      LO_Cut,
      LO_Link,
      LO_Vertex,
      LO_Cross,
     )

type LocalOperator
    op_type :: LocalOperatorType
    isdiagonal :: Bool
    time :: Float64
    space :: Int
    bottom_id :: Int
    top_id :: Int
end
LocalOperator(op_type::LocalOperatorType, time::Real, space::Int) = LocalOperator(op_type, true, time, space, 0,0)

abstract QuantumLocalZ2Model <: Model

type TransverseFieldIsing <: QuantumLocalZ2Model
    lat :: Lattice
    spins :: Vector{Int}
    ops :: Vector{LocalOperator}

    function TransverseFieldIsing(lat::Lattice)
        model = new()
        model.lat = lat
        model.spins = rand([1,-1], numsites(lat))
        model.ops = LocalOperator[]
        return model
    end
end
function TransverseFieldIsing(params::Dict)
    lat = params["Lattice"](params)
    return TransverseFieldIsing(lat)
end

