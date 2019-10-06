__precompile__()

module SpinMonteCarlo

using Random
using Printf
using Markdown
using Statistics
using LinearAlgebra

using DataStructures

@doc doc"""
Input parameter of simulation
"""
const Measurement = Dict{String, Any}

include("types.jl")

include("lattice/Lattices.jl")
include("model/model.jl")
include("union_find.jl")
include("update/update.jl")
include("estimator/estimator.jl")
include("snapshot.jl")
include("observables/MCObservables.jl")
include("parameter.jl")
include("runMC.jl")
include("print.jl")
end # of module
