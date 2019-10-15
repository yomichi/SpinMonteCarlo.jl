__precompile__()

module SpinMonteCarlo

using Random
using Printf
using Markdown
using Statistics
using LinearAlgebra

using DataStructures

include("observables/MCObservables.jl")
include("API/api.jl")
include("model/model.jl")
include("lattice/Lattices.jl")
include("runMC.jl")
include("snapshot.jl")
end # of module
