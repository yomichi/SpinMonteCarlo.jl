__precompile__()

module SpinMonteCarlo

using Random
using Printf
using Markdown
using Statistics
using LinearAlgebra

using DataStructures

include("API/api.jl")
include("model/model.jl")
include("lattice/Lattices.jl")
include("observables/MCObservables.jl")
include("runMC.jl")
include("snapshot.jl")
end # of module
