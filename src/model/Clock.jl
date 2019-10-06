@doc doc"""
`Q` state clock model with energy $E = -\sum_{ij} J_{ij} \cos(\theta_i - \theta_j)$,
where $\theta_i = 2\pi \sigma_i/Q$ and $\sigma_i$ takes an integer value from $1$ to $Q$.
"""
mutable struct Clock <: Model
    lat :: Lattice
    Q :: Int
    spins :: Matrix{Int}
    cosines :: Vector{Float64}
    sines :: Vector{Float64}
    sines_sw :: Vector{Float64}
    rng :: Random.MersenneTwister

    function Clock(lat::Lattice, Q::Integer)
        rng = Random.seed!(Random.MersenneTwister(0))
        spins = rand(rng, 1:Q, 1, numsites(lat))
        cosines = [cospi(2s/Q) for s in 1:Q]
        sines = [sinpi(2s/Q) for s in 1:Q]
        sines_sw = [sinpi(2(s-0.5)/Q) for s in 1:Q]
        return new(lat, Q, spins, cosines, sines, sines_sw, rng)
    end
    function Clock(lat::Lattice, Q::Integer, seed)
        rng = Random.seed!(Random.MersenneTwister(0), seed)
        spins = rand(rng, 1:Q, 1, numsites(lat))
        cosines = [cospi(2s/Q) for s in 1:Q]
        sines = [sinpi(2s/Q) for s in 1:Q]
        sines_sw = [sinpi(2(s-0.5)/Q) for s in 1:Q]
        return new(lat, Q, spins, cosines, sines, sines_sw, rng)
    end
end
@doc doc"""
    Clock(param)

Generates `Clock` using `param["Lattice"]`, `param["Q"]`,  and `param["Seed"]` (if defined).
Each spin $\sigma_i$ will be initialized randomly and independently.
"""
function Clock(param::Parameter)
    lat = generatelattice(param)
    Q = param["Q"]
    if "Seed" in keys(param)
        return Clock(lat, Q, param["Seed"])
    else
        return Clock(lat, Q)
    end
end
