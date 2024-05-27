__precompile__()

module SpinMonteCarlo

@warn "From 2024-05-27, the default branch of SpinMonteCarlo.jl is 'main', but you are in the old branch 'master'. Please migrate to 'main'."

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
