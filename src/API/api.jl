export Model, Parameter, Measurement

abstract type Model end


@doc doc"""
Input parameter of simulation
"""
const Parameter = Dict{String,Any}

const Measurement = Dict{String, Any}

include("lattice.jl")
include("update.jl")
include("estimator.jl")
include("parameter.jl")
include("postproc.jl")
