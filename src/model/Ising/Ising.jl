@doc doc"""
Ising model with energy $E = -\sum_{ij} J_{ij} \sigma_i \sigma_j$,
where $\sigma_i$ takes value of 1 (up spin) or -1 (down spin).
"""
mutable struct Ising <: Model
    lat :: Lattice
    spins :: Matrix{Int}
    rng :: Random.MersenneTwister

    function Ising(lat::Lattice, rng::Random.AbstractRNG)
        model = new()
        model.lat = lat
        model.rng = rng
        model.spins = rand(model.rng, [1,-1], 1, numsites(lat))
        return model
    end
end

Ising(lat::Lattice) = Ising(lat, Random.seed!(Random.MersenneTwister(0)))
Ising(lat::Lattice, seed) = Ising(lat, Random.seed!(Random.MersenneTwister(0), seed...))

@doc doc"""
    Ising(param)

Generates `Ising` using `param["Lattice"]` and `param["Seed"]` (if defined).
Each spin will be initialized randomly and independently.
"""
function Ising(param::Parameter)
    lat = generatelattice(param)
    if "Seed" in keys(param)
        return Ising(lat, param["Seed"])
    else
        return Ising(lat)
    end
end

include("update.jl")
include("estimator.jl")
include("postproc.jl")
