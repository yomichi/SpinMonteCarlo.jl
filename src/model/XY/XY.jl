@doc doc"""
XY model with energy $E = -\sum_{ij} J_{ij} \cos(\theta_i - \theta_j)$,
where $\theta_i = 2\pi \sigma_i$ and $\sigma_i \in [0, 1)$.
"""
mutable struct XY{RNG<:Random.AbstractRNG} <: Model
    lat::Lattice
    spins::Matrix{Float64}
    rng::RNG

    function XY(lat::Lattice, rng::R) where {R<:Random.AbstractRNG}
        spins = rand(rng, 1, numsites(lat))
        return new{R}(lat, spins, rng)
    end
end

XY(lat::Lattice) = XY(lat, DEFAULT_RNG())
XY(lat::Lattice, seed) = XY(lat, DEFAULT_RNG(seed))

@doc doc"""
   XY(param)

Generates `XY` using `param["Lattice"]`, and `param["Seed"]` and `param["RNG"]` (if defined).
Each spin $\sigma_i$ will be initialized randomly and independently.
"""
function XY(param::Parameter)
    lat = generatelattice(param)
    return XY(lat, makerng(param))
end

include("update.jl")
include("estimator.jl")
include("postproc.jl")
