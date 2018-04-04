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

type QuantumXXZ <: QuantumLocalZ2Model
    lat :: Lattice
    S2 :: Int
    spins :: Vector{Int}
    ops :: Vector{LocalOperator}

    function QuantumXXZ(lat::Lattice, S2::Int)
        model = new()
        model.lat = lat
        model.S2 = S2
        model.spins = rand([1,-1], numsites(lat)*S2)
        model.ops = LocalOperator[]
        return model
    end
end
function QuantumXXZ(params::Dict)
    lat = params["Lattice"](params)
    S = params["S"]
    if round(2S) != 2S
        error("`S` should be integer or half-integer")
    end
    S2 = round(Int,2S)
    return QuantumXXZ(lat,S2)
end

site2subspin(site::Integer, ss::Integer, S2::Integer) = (site-1)*S2+ss
bond2subbond(bond::Integer, ss1::Integer, ss2::Integer, S2::Integer) = ((bond-1)*S2+(ss1-1))*S2+ss2
function subspin2site(subspin::Integer, S2::Integer)
    return ceil(Int, subspin/S2), mod1(subspin,S2)
end
function subbond2bond(subbond::Integer, S2::Integer)
    ss2 = mod1(subbond,S2)
    ss = ceil(Int, subbond/S2)
    return ceil(Int, ss/S2), mod1(ss,S2), ss2
end
