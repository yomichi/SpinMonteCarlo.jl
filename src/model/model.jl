import Random.seed!
seed!(model::Model) = Random.seed!(model.rng)
seed!(model::Model, seed...) = Random.seed!(model.rng, seed...)

export Model, Ising, XY, Potts, Clock
export QuantumLocalZ2Model, QuantumXXZ

include("Ising.jl")
include("Potts.jl")
include("Clock.jl")
include("XY.jl")
include("qmodel.jl")
