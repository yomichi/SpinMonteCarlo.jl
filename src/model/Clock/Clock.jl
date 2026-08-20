@doc doc"""
`Q` state clock model with energy $E = -\sum_{ij} J_{ij} \cos(\theta_i - \theta_j)$,
where $\theta_i = 2\pi \sigma_i/Q$ and $\sigma_i$ takes an integer value from $1$ to $Q$.
"""
mutable struct Clock{RNG<:Random.AbstractRNG} <: Model
    lat::Lattice
    Q::Int
    spins::Matrix{Int}
    cosines::Vector{Float64}
    sines::Vector{Float64}
    sines_sw::Vector{Float64}
    rng::RNG

    function Clock(lat::Lattice, Q::Integer, rng::R) where {R<:Random.AbstractRNG}
        if Q < 2
            error("Q should be 2 or more: Q = $Q")
        end
        spins = rand(rng, 1:Q, 1, numsites(lat))
        cosines = [cospi(2s / Q) for s in 1:Q]
        sines = [sinpi(2s / Q) for s in 1:Q]
        sines_sw = [sinpi(2(s - 0.5) / Q) for s in 1:Q]
        return new{R}(lat, Q, spins, cosines, sines, sines_sw, rng)
    end
end
Clock(lat::Lattice, Q::Integer) = Clock(lat, Q, DEFAULT_RNG())
function Clock(lat::Lattice, Q::Integer, seed)
    return Clock(lat, Q, DEFAULT_RNG(seed))
end

@doc doc"""
    Clock(param)

Generates `Clock` using `param["Lattice"]` and `param["Q"]`, and `param["Seed"]` and `param["RNG"]` (if defined).
Each spin $\sigma_i$ will be initialized randomly and independently.
"""
function Clock(param::Parameter)
    lat = generatelattice(param)
    Q = param["Q"]
    return Clock(lat, Q, makerng(param))
end

include("update.jl")
include("estimator.jl")
include("postproc.jl")
