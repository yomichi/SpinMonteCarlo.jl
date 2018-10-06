using LightGraphs
using MetaGraphs

const Lattice = MetaGraph{Int,Float64}

include("parameter.jl")
include("api.jl")
include("generate.jl")
include("standard.jl")

export LatticeParameter
export generatelattice
