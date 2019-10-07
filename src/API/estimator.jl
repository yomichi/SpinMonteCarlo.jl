export simple_estimator, improved_estimator, default_estimator

@inline function simple_estimator(model, param::Parameter, _=nothing)
    p = convert_parameter(model,param)
    return simple_estimator(model, p..., nothing)
end

@inline function improved_estimator(model, param::Parameter, extra)
    p = convert_parameter(model, param)
    return improved_estimator(model, p..., extra)
end

@doc doc"""
    default_estimator(model, updatemethod!)

Determines estimator to be used when `param["Estimator"]` is not set.
"""
@inline default_estimator(model::Model, update) = simple_estimator
