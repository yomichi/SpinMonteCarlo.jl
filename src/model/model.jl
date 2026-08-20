import Random.seed!
seed!(model::Model) = Random.seed!(model.rng)
seed!(model::Model, seed) = Random.seed!(model.rng, seed)

const DEFAULT_RNG = Random.Xoshiro
const SEED_GAMMA = 0x9e3779b97f4a7c15

function splitmix64_finalize(x::UInt64)
    x ⊻= x >> 30
    x *= 0xbf58476d1ce4e5b9
    x ⊻= x >> 27
    x *= 0x94d049bb133111eb
    x ⊻= x >> 31
    return x
end

@doc doc"""
    childseed(seed::Integer, id::Integer)

Derive a deterministic child seed from an integer seed and an integer ID.

This uses a fixed SplitMix64 finalizer rather than `Base.hash`, so integer
seed derivation is stable across Julia versions. Child seed derivation is
used by `makerng(param)` only when `param` has an `"ID"` key, such as IDs set
by `runMC(params)` with `autoID=true`. Passing a model directly to
`runMC(model, param)` does not derive a child seed.
"""
function childseed(seed::Integer, id::Integer)
    return splitmix64_finalize(splitmix64_finalize((seed % UInt64) + SEED_GAMMA) ⊻
                               (SEED_GAMMA * (id % UInt64)))
end

seedvalue(seed::Integer) = seed
seedvalue(seed) = hash(seed)

@doc doc"""
    makerng(param::Parameter)

Create the random number generator requested by `param`.

`param["RNG"]`, when present, must be an RNG type, not an RNG instance. It
needs only the constructor the other keys call for: `T(seed)` when
`param["Seed"]` is given, and `T()` when it is not. `Random.RandomDevice` and
`Random.TaskLocalRNG`, for instance, provide the latter but not the former, so
they work only without `"Seed"`.

When `param` has both `"Seed"` and `"ID"`, a child seed is derived for that
ID. `runMC(params)` sets these IDs when `autoID=true`; `runMC(model, param)`
does not derive child seeds because the model already owns its RNG. Integer
seeds use a fixed mixer and are reproducible across Julia versions. Non-integer
seeds combined with `"ID"` pass through `hash`, so only that case is not
guaranteed to be reproducible across Julia versions.
"""
function makerng(param::Parameter)
    RNG = get(param, "RNG", DEFAULT_RNG)
    haskey(param, "Seed") || return RNG()
    seed = param["Seed"]
    haskey(param, "ID") || return RNG(seed)
    return RNG(childseed(seedvalue(seed), param["ID"]))
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
