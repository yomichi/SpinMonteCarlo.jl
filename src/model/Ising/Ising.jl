@doc doc"""
Ising model with energy $E = -\sum_{ij} J_{ij} \sigma_i \sigma_j$,
where $\sigma_i$ takes value of 1 (up spin) or -1 (down spin).
"""
mutable struct Ising{RNG<:Random.AbstractRNG} <: Model
    lat::Lattice
    spins::Matrix{Int}
    rng::RNG

    function Ising(lat::Lattice, rng::R) where {R<:Random.AbstractRNG}
        model = new{R}()
        model.lat = lat
        model.rng = rng
        model.spins = rand(model.rng, [1, -1], 1, numsites(lat))
        return model
    end
end

Ising(lat::Lattice) = Ising(lat, DEFAULT_RNG())
Ising(lat::Lattice, seed) = Ising(lat, DEFAULT_RNG(seed))

@doc doc"""
    Ising(param)

Generates `Ising` using `param["Lattice"]`, and `param["Seed"]` and `param["RNG"]` (if defined).
Each spin will be initialized randomly and independently.
"""
function Ising(param::Parameter)
    lat = generatelattice(param)
    return Ising(lat, makerng(param))
end

include("update.jl")
include("estimator.jl")
include("postproc.jl")
