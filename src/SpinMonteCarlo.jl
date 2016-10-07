module SpinMonteCarlo

export Ising, XY, Model
export SW_update!, local_update!
export chain_lattice, square_lattice, triangular_lattice
export Lattice, dim, size, num_sites, num_bonds, neighbors, source, target

include("union_find.jl")
include("lattice.jl")
include("model.jl")
include("SW.jl")
include("local_update.jl")
include("observables/MCObservables.jl")

end # of module
