@doc doc"""
`Q` state Potts model with energy $E = -\sum_{i,j} \delta_{\sigma_i, \sigma_j}$,
where $\sigma_i$ takes an integer value from $1$ to $Q$ and $\delta$ is a Kronecker's delta.
Order parameter (total magnetization) is defined as
\begin{equation}
    M = \frac{Q-1}{Q}N_1 - \frac{1}{Q}(N-N_1),
\end{equation}
where $N$ is the number of sites and $N_1$ is the number of $\sigma=1$ spins.
"""
mutable struct Potts{RNG<:Random.AbstractRNG} <: Model
    lat::Lattice
    Q::Int
    spins::Matrix{Int}
    rng::RNG

    function Potts(lat::Lattice, Q::Integer, rng::R) where {R<:Random.AbstractRNG}
        if Q < 2
            error("Q should be 2 or more: Q = $Q")
        end
        spins = rand(rng, 1:Q, 1, numsites(lat))
        return new{R}(lat, Q, spins, rng)
    end
end

Potts(lat::Lattice, Q::Integer) = Potts(lat, Q, DEFAULT_RNG())
function Potts(lat::Lattice, Q::Integer, seed)
    return Potts(lat, Q, DEFAULT_RNG(seed))
end

@doc doc"""
    Potts(param)

Generates `Potts` using `param["Lattice"]`, `param["Q"]`, and `param["Seed"]` (if defined).
Each spin will be initialized randomly and independently.
"""
function Potts(param::Parameter)
    lat = generatelattice(param)
    Q = param["Q"]
    return Potts(lat, Q, makerng(param))
end

include("update.jl")
include("estimator.jl")
include("postproc.jl")
