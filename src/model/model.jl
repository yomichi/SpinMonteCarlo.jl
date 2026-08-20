import Random.seed!
seed!(model::Model) = Random.seed!(model.rng)
seed!(model::Model, seed...) = Random.seed!(model.rng, seed...)

const DEFAULT_RNG = Random.Xoshiro
@doc doc"""
    makerng(param::Parameter)

Create the random number generator requested by `param`.

`param["RNG"]`, when present, must be an RNG type, not an RNG instance. If
`param["Seed"]` is also present, that type must support both `T()` and
`T(seed)` constructors. For example, `Random.RandomDevice` and
`Random.TaskLocalRNG` cannot be used together with `"Seed"`.
"""
function makerng(param::Parameter)
    RNG = get(param, "RNG", DEFAULT_RNG)
    haskey(param, "Seed") || return RNG()
    return RNG(param["Seed"])
end

export Model, Ising, XY, Potts, Clock, AshkinTeller
export QuantumXXZ

include("common/union_find.jl")

## Classical
include("Ising/Ising.jl")
include("Potts/Potts.jl")
include("Clock/Clock.jl")
include("XY/XY.jl")
include("AshkinTeller/AshkinTeller.jl")

## Quantum
include("common/LoopOperator.jl")
include("QuantumXXZ/QuantumXXZ.jl")
