export default_estimator, simple_estimator, improved_estimator

@doc doc"""
    default_estimator(model, updatemethod!)

Determines estimator to be used when `param["Estimator"]` is not set.
"""
default_estimator(model::Model, update) = simple_estimator
default_estimator(model::QuantumXXZ, update) = improved_estimator
default_estimator(model::Union{Ising, Potts}, update) = ifelse(update==SW_update!, improved_estimator, simple_estimator)


include("simple_estimator.jl")
include("improved_estimator.jl")
