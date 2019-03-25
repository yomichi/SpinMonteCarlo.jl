import Random.seed!
seed!(model::Model) = Random.seed!(model.rng)
seed!(model::Model, seed...) = Random.seed!(model.rng, seed...)

@doc doc"""
Ising model with energy $E = -\sum_{ij} J_{ij} \sigma_i \sigma_j$,
where $\sigma_i$ takes value of 1 (up spin) or -1 (down spin).
"""
mutable struct Ising <: Model
    lat :: Lattice
    spins :: Vector{Int}
    rng :: Random.MersenneTwister
    parameter :: Parameter

    function Ising(lat::Lattice)
        model = new()
        model.lat = lat
        model.rng = Random.Random.seed!(Random.MersenneTwister(0))
        model.spins = rand(model.rng, [1,-1], numsites(lat))
        model.parameter = Parameter()
        return model
    end
    function Ising(lat::Lattice, seed)
        model = new()
        model.lat = lat
        model.rng = Random.seed!(Random.MersenneTwister(0), seed...)
        model.spins = rand(model.rng, [1,-1], numsites(lat))
        model.parameter = Parameter()
        return model
    end
end
@doc doc"""
    Ising(param)

Generates `Ising` using `param["Lattice"]` and `param["Seed"]` (if defined).
Each spin will be initialized randomly and independently.
"""
function Ising(param::Parameter, lat::Lattice)
    if "Seed" in keys(param)
        model = Ising(lat, param["Seed"])
    else
        model = Ising(lat)
    end
    merge!(model.parameter, param)
    return model
end
Ising(param::Parameter) = Ising(param, generatelattice(param))

@doc doc"""
`Q` state Potts model with energy $E = -\sum_{i,j} \delta_{\sigma_i, \sigma_j}$,
where $\sigma_i$ takes an integer value from $1$ to $Q$ and $\delta$ is a Kronecker's delta.
Order parameter (total magnetization) is defined as
\begin{equation}
    M = \frac{Q-1}{Q}N_1 - \frac{1}{Q}(N-N_1),
\end{equation}
where $N$ is the number of sites and $N_1$ is the number of $\sigma=1$ spins.
"""
mutable struct Potts <: Model
    lat :: Lattice
    Q :: Int
    spins :: Vector{Int}
    rng :: Random.MersenneTwister
    parameter :: Parameter

    function Potts(lat::Lattice, Q::Integer)
        rng = Random.seed!(Random.MersenneTwister(0))
        spins = rand(rng, 1:Q, numsites(lat))
        return new(lat, Q, spins, rng, Parameter())
    end
    function Potts(lat::Lattice, Q::Integer, seed)
        rng = Random.seed!(Random.MersenneTwister(0), seed)
        spins = rand(rng, 1:Q, numsites(lat))
        return new(lat, Q, spins, rng, Parameter())
    end
end
@doc doc"""
    Potts(param)

Generates `Potts` using `param["Lattice"]`, `param["Q"]`, and `param["Seed"]` (if defined).
Each spin will be initialized randomly and independently.
"""
function Potts(param::Parameter, lat::Lattice)
    Q = param["Q"]
    if "Seed" in keys(param)
        model = Potts(lat, Q, param["Seed"])
    else
        model = Potts(lat, Q)
    end
    merge!(model.parameter, param)
    return model
end
Potts(param::Parameter) = Potts(param, generatelattice(param))

@doc doc"""
`Q` state clock model with energy $E = -\sum_{ij} J_{ij} \cos(\theta_i - \theta_j)$,
where $\theta_i = 2\pi \sigma_i/Q$ and $\sigma_i$ takes an integer value from $1$ to $Q$.
"""
mutable struct Clock <: Model
    lat :: Lattice
    Q :: Int
    spins :: Vector{Int}
    cosines :: Vector{Float64}
    sines :: Vector{Float64}
    sines_sw :: Vector{Float64}
    rng :: Random.MersenneTwister
    parameter :: Parameter

    function Clock(lat::Lattice, Q::Integer)
        rng = Random.seed!(Random.MersenneTwister(0))
        spins = rand(rng, 1:Q, numsites(lat))
        cosines = [cospi(2s/Q) for s in 1:Q]
        sines = [sinpi(2s/Q) for s in 1:Q]
        sines_sw = [sinpi(2(s-0.5)/Q) for s in 1:Q]
        return new(lat, Q, spins, cosines, sines, sines_sw, rng, Parameter())
    end
    function Clock(lat::Lattice, Q::Integer, seed)
        rng = Random.seed!(Random.MersenneTwister(0), seed)
        spins = rand(rng, 1:Q, numsites(lat))
        cosines = [cospi(2s/Q) for s in 1:Q]
        sines = [sinpi(2s/Q) for s in 1:Q]
        sines_sw = [sinpi(2(s-0.5)/Q) for s in 1:Q]
        return new(lat, Q, spins, cosines, sines, sines_sw, rng, Parameter())
    end
end
@doc doc"""
    Clock(param)

Generates `Clock` using `param["Lattice"]`, `param["Q"]`,  and `param["Seed"]` (if defined).
Each spin $\sigma_i$ will be initialized randomly and independently.
"""
function Clock(param::Parameter, lat::Lattice)
    Q = param["Q"]
    if "Seed" in keys(param)
        model = Clock(lat, Q, param["Seed"])
    else
        model = Clock(lat, Q)
    end
    merge!(model.parameter, param)
    return model
end
Clock(param::Parameter) = Clock(param, generatelattice(param))

@doc doc"""
XY model with energy $E = -\sum_{ij} J_{ij} \cos(\theta_i - \theta_j)$,
where $\theta_i = 2\pi \sigma_i$ and $\sigma_i \in [0, 1)$.
"""
mutable struct XY <: Model
    lat :: Lattice
    spins :: Vector{Float64}
    rng :: Random.MersenneTwister
    parameter :: Parameter

    function XY(lat::Lattice)
        model = new()
        model.rng = Random.seed!(Random.MersenneTwister(0))
        model.lat = lat
        model.spins = rand(model.rng, numsites(lat))
        model.parameter = Parameter()
        return model
    end
    function XY(lat::Lattice, seed)
        model = new()
        model.rng = Random.seed!(Random.MersenneTwister(0), seed)
        model.lat = lat
        model.spins = rand(model.rng, numsites(lat))
        model.parameter = Parameter()
        return model
    end
end
@doc doc"""
   XY(param)

Generates `XY` using `param["Lattice"]`,  and `param["Seed"]` (if defined).
Each spin $\sigma_i$ will be initialized randomly and independently.
"""
function XY(param::Parameter, lat::Lattice)
    if "Seed" in keys(param)
        model = XY(lat, param["Seed"])
    else
        model = XY(lat)
    end
    merge!(model.parameter, param)
    return model
end
XY(param::Parameter) = XY(param, generatelattice(param))

