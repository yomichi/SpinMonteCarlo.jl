@doc doc"""
XY model with energy $E = -\sum_{ij} J_{ij} \cos(\theta_i - \theta_j)$,
where $\theta_i = 2\pi \sigma_i$ and $\sigma_i \in [0, 1)$.
"""
mutable struct XY <: Model
    lat :: Lattice
    spins :: Matrix{Float64}
    rng :: Random.MersenneTwister

    function XY(lat::Lattice)
        model = new()
        model.rng = Random.seed!(Random.MersenneTwister(0))
        model.lat = lat
        model.spins = rand(model.rng, 1, numsites(lat))
        return model
    end
    function XY(lat::Lattice, seed)
        model = new()
        model.rng = Random.seed!(Random.MersenneTwister(0), seed)
        model.lat = lat
        model.spins = rand(model.rng, 1, numsites(lat))
        return model
    end
end
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
