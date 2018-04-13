@enum(LocalOperatorType,
      LO_Cut,    # [1 1; 1 1]
      LO_FMLink, # [1 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 1]
      LO_AFLink, # [0 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]
      LO_Vertex, # [0 0 0 0; 0 1 1 0; 0 1 1 0; 0 0 0 0]
      LO_Cross,  # [1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1]
     )

mutable struct LocalOperator
    op_type :: LocalOperatorType
    isdiagonal :: Bool
    time :: Float64
    space :: Int
    bottom_id :: Int
    top_id :: Int
end
LocalOperator(op_type::LocalOperatorType, time::Real, space::Int) = LocalOperator(op_type, true, time, space, 0,0)

abstract type QuantumLocalZ2Model <: Model end

doc"""
    QuantumXXZ(lat::Lattice, S::Real, [seed])

Spin-$S$ XXZ model denoted by the following Hamiltonian,
\begin{equation}
\mathcal{H} = \sum_{i,j} \left[ J_{ij}^z S_i^z S_j^z 
            + \frac{J_{ij}^{xy}}{2} (S_i^+ S_j^- + S_i^-S_j^+) \right]
            - \sum_i \Gamma_i S_i^x,
\end{equation}
where $S^x, S^y, S^z$ are $x, y$ and $z$ component of spin operator with length $S$,
and $S^\pm \equiv S^x \pm iS^y$ are ladder operator.
A state is represented by a product state (spins at $\tau=0$) of local $S^z$ diagonal basis and an operator string (perturbations).
A local spin with length $S$ is represented by a symmetrical summation of $2S$ sub spins with length $1/2$.
Each subspin will be initialized independently and randomly.

    QuantumXXZ(param)

`param["Lattice"](param)`, `param["S"]` and `param["Seed"]` (if defined) will be used as the args of the former form.
"""
mutable struct QuantumXXZ <: QuantumLocalZ2Model
    lat :: Lattice
    S2 :: Int
    spins :: Vector{Int}
    ops :: Vector{LocalOperator}
    rng :: Random.MersenneTwister

    function QuantumXXZ(lat::Lattice, S::Real)
        if round(2S) != 2S
            error("`S` should be integer or half-integer")
        end
        S2 = round(Int,2S)
        model = new()
        model.rng = srand(Random.MersenneTwister(0))
        model.lat = lat
        model.S2 = S2
        model.spins = rand(model.rng,[1,-1], numsites(lat)*S2)
        model.ops = LocalOperator[]
        return model
    end
    function QuantumXXZ(lat::Lattice, S::Real, seed)
        if round(2S) != 2S
            error("`S` should be integer or half-integer")
        end
        S2 = round(Int,2S)
        model = new()
        model.rng = srand(Random.MersenneTwister(0), seed)
        model.lat = lat
        model.S2 = S2
        model.spins = rand(model.rng,[1,-1], numsites(lat)*S2)
        model.ops = LocalOperator[]
        return model
    end
end
function QuantumXXZ(param::Dict)
    lat = param["Lattice"](param)
    S = param["S"]
    if "Seed" in keys(param)
        return QuantumXXZ(lat,S,param["Seed"])
    else
        return QuantumXXZ(lat,S)
    end
end

site2subspin(site::Integer, ss::Integer, S2::Integer) = (site-1)*S2+ss
bond2subbond(bond::Integer, ss1::Integer, ss2::Integer, S2::Integer) = ((bond-1)*S2+(ss1-1))*S2+ss2
@inline function subspin2site(subspin::Integer, S2::Integer)
    return ceil(Int, subspin/S2), mod1(subspin,S2)
end
@inline function subbond2bond(subbond::Integer, S2::Integer)
    ss2 = mod1(subbond,S2)
    ss = ceil(Int, subbond/S2)
    return ceil(Int, ss/S2), mod1(ss,S2), ss2
end
