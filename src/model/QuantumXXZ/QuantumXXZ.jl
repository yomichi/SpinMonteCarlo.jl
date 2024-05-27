@doc raw"""
Spin-``S`` XXZ model represented as the following Hamiltonian,

``
\mathcal{H} = \sum_{i,j} \left[ J_{ij}^z S_i^z S_j^z 
            + \frac{J_{ij}^{xy}}{2} (S_i^+ S_j^- + S_i^-S_j^+) \right]
            - \sum_i \Gamma_i S_i^x,
``

where ``S^x, S^y, S^z`` are ``x, y`` and ``z`` component of spin operator with length ``S``,
and ``S^\pm \equiv S^x \pm iS^y`` are ladder operator.
A state is represented by a product state (spins at ``\tau=0``) of local ``S^z`` diagonal basis and an operator string (perturbations).
A local spin with length ``S`` is represented by a symmetrical summation of ``2S`` sub spins with length ``1/2``.
"""
mutable struct QuantumXXZ <: Model
    lat::Lattice
    S2::Int
    spins::Matrix{Int}
    ops::Vector{LocalLoopOperator}
    rng::Random.MersenneTwister

    function QuantumXXZ(lat::Lattice, S::Real, rng::Random.AbstractRNG)
        if round(2S) != 2S
            error("`S` should be integer or half-integer")
        end
        S2 = round(Int, 2S)
        model = new()
        model.lat = lat
        model.rng = rng
        model.S2 = S2
        model.spins = rand(model.rng, [1, -1], 1, numsites(lat) * S2)
        model.ops = LocalLoopOperator[]
        return model
    end
end
function QuantumXXZ(lat::Lattice, S::Real)
    return QuantumXXZ(lat, S, Random.seed!(Random.MersenneTwister(0)))
end
function QuantumXXZ(lat::Lattice, S::Real, seed)
    return QuantumXXZ(lat, S, Random.seed!(Random.MersenneTwister(0), seed...))
end

@doc doc"""
    QuantumXXZ(param)

Generates `QuantumXXZ` using `param["Lattice"]`, `param["S"]` and `param["Seed"]` (if defined).
Each subspin will be initialized independently and randomly.
"""
function QuantumXXZ(param::Parameter)
    lat = generatelattice(param)
    S = param["S"]
    if "Seed" in keys(param)
        return QuantumXXZ(lat, S, param["Seed"])
    else
        return QuantumXXZ(lat, S)
    end
end

@inline site2subspin(site::Integer, ss::Integer, S2::Integer) = (site - 1) * S2 + ss
@inline function subspin2site(subspin::Integer, S2::Integer)
    return ceil(Int, subspin / S2), mod1(subspin, S2)
end

include("update.jl")
include("estimator.jl")
include("postproc.jl")
