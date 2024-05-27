@doc doc"""
XY model with energy $E = -\sum_{ij} J_{ij} \cos(\theta_i - \theta_j)$,
where $\theta_i = 2\pi \sigma_i$ and $\sigma_i \in [0, 1)$.
"""
mutable struct XY <: Model
    lat::Lattice
    spins::Matrix{Float64}
    rng::Random.MersenneTwister

    function XY(lat::Lattice, rng::Random.AbstractRNG)
        spins = rand(rng, 1, numsites(lat))
        return new(lat, spins, rng)
    end
end

XY(lat::Lattice) = XY(lat, Random.seed!(Random.MersenneTwister(0)))
XY(lat::Lattice, seed) = XY(lat, Random.seed!(Random.MersenneTwister(0), seed...))

@doc doc"""
   XY(param)

Generates `XY` using `param["Lattice"]`,  and `param["Seed"]` (if defined).
Each spin $\sigma_i$ will be initialized randomly and independently.
"""
function XY(param::Parameter)
    lat = generatelattice(param)
    if "Seed" in keys(param)
        return XY(lat, param["Seed"])
    else
        return XY(lat)
    end
end

include("update.jl")
include("estimator.jl")
include("postproc.jl")
