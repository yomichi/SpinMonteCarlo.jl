@doc raw"""
AshkinTeller model with energy
``E = -\sum_{ij} J^\sigma_{ij} \sigma_i \sigma_j + J^\tau_{ij} \tau_i \tau_j + K_{ij} \sigma_i \sigma_j \tau_i \tau_j``,
where ``\sigma_i`` and ``\tau_i`` takes value of 1 (up spin) or -1 (down spin).
"""
mutable struct AshkinTeller <: Model
    lat::Lattice
    spins::Matrix{Int}
    rng::Random.MersenneTwister

    function AshkinTeller(lat::Lattice, rng::Random.AbstractRNG)
        model = new()
        model.lat = lat
        model.rng = rng
        model.spins = rand(model.rng, [1, -1], 2, numsites(lat))
        return model
    end
end

AshkinTeller(lat::Lattice) = AshkinTeller(lat, Random.seed!(Random.MersenneTwister(0)))
function AshkinTeller(lat::Lattice, seed)
    return AshkinTeller(lat, Random.seed!(Random.MersenneTwister(0), seed...))
end

@doc doc"""
    AshkinTeller(param)

Generates `AshkinTeller` using `param["Lattice"]` and `param["Seed"]` (if defined).
Each spin will be initialized randomly and independently.
"""
function AshkinTeller(param::Parameter)
    lat = generatelattice(param)
    if "Seed" in keys(param)
        return AshkinTeller(lat, param["Seed"])
    else
        return AshkinTeller(lat)
    end
end

include("update.jl")
include("estimator.jl")
include("postproc.jl")
