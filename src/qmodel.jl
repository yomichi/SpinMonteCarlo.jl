@enum(LocalOperatorType,
      LO_Cut,    # [1 1; 1 1]
      LO_FMLink, # [1 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 1]
      LO_AFLink, # [0 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]
      LO_Vertex, # [0 0 0 0; 0 1 1 0; 0 1 1 0; 0 0 0 0]
      LO_Cross,  # [1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1]
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

@compat abstract type QuantumLocalZ2Model <: Model end

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

