VERSION >= v"0.4.0-dev+6521" && __precompile__()

module SpinMonteCarlo
using Compat
using DataStructures

export Model, Ising, XY, Potts, Clock
export QuantumLocalZ2Model, TransverseFieldIsing
export local_update!, SW_update!, Wolff_update!, loop_update!
export magnetizations, energy
export chain_lattice, square_lattice, triangular_lattice, cubic_lattice
export Lattice, dim, size, numsites, numbonds, neighbors, source, target, siteL2, siteL4
export UnionFind, addnode!, unify!, clusterize!, clusterid
export measure
export gen_snapshot!, gensave_snapshot!, load_snapshot
export runMC, print_result

const Measurement = Dict{Symbol, Any}

include("union_find.jl")
include("lattice.jl")
include("model.jl")
include("qmodel.jl")
include("SW.jl")
include("loop.jl")
include("local_update.jl")
include("Wolff.jl")
include("simple_estimator.jl")
include("improved_estimator.jl")
include("snapshot.jl")
include("observables/MCObservables.jl")
include("runMC.jl")
include("print.jl")

end # of module
