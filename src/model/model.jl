import Random.seed!
seed!(model::Model) = Random.seed!(model.rng)
seed!(model::Model, seed...) = Random.seed!(model.rng, seed...)

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
