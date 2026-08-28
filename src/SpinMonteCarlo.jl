module SpinMonteCarlo

using Random
using Printf
using Markdown
using Statistics
using LinearAlgebra

include("observables/MCObservables.jl")
include("API/api.jl")
include("model/model.jl")
include("lattice/Lattices.jl")
include("snapshot.jl")
include("runMC.jl")
end # of module
