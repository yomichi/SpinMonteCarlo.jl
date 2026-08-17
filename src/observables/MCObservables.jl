using SpecialFunctions
using LsqFit

import Distributions: TDist, cdf

import Base: show, <<, push!, count, isempty, merge, merge!, zero, zeros, sum

import Statistics: mean, var

export MCObservable, ScalarObservable, VectorObservable
export mean, var, stderror, confidence_interval
export p_value
export show, dump_plot
export merge, merge!

@doc """
    MCObservable

    abstract type representing observable in Monte Carlo calculation.
"""
abstract type MCObservable end

"""
    ScalarObservable <: MCObservable

    abstract type representing scalar-type observable in Monte Carlo calculation.
"""
abstract type ScalarObservable <: MCObservable end

"""
    VectorObservable <: MCObservable

    abstract type representing vector-type observable in Monte Carlo calculation.
"""
abstract type VectorObservable <: MCObservable end

function p_value(X::MCObservable, y::Real; verbose::Bool=false)
    MX = mean(X)
    mid = 0.5(MX + y)
    er = max(stderror(X), abs(mid) * sqrt(eps()), sqrt(nextfloat(0.0)))
    t = (MX - y) / er
    d = TDist(count(X) - 1)
    if verbose
        @show MX, y, er, t
    end
    return 2cdf(d, -abs(t))
end

function p_value(X::MCObservable, Y::MCObservable)
    NX = count(X)
    NY = count(Y)
    MX = mean(X)
    MY = mean(Y)
    mid = 0.5(MX + MY)
    N = 1.0 / (1.0 / NX + 1.0 / NY)
    V = ((NX - 1) * var(X) + (NY - 1) * var(Y)) / (NX + NY - 2)
    se = max(sqrt(V / N), abs(mid) * sqrt(eps()), sqrt(nextfloat(0.0)))
    t = (MX - MY) / se
    d = TDist(NX + NY - 2)
    return 2cdf(d, -abs(t))
end

isempty(obs::MCObservable) = count(obs) == 0

zero(::Type{Obs}) where {(Obs<:MCObservable)} = Obs()
zero(obs::MCObservable) = zero(typeof(obs))

## Dimensions are taken as a tuple rather than as a vararg, since a vararg
## method is ambiguous with the `zeros(::Type, ::Tuple)` methods of Base.
## Note that each element has to be a distinct object, so `fill`, which Base
## uses, is not applicable here.
function zeros(::Type{Obs}, dims::Tuple{Vararg{Integer}}) where {(Obs<:MCObservable)}
    return reshape(Obs[Obs() for i in 1:prod(dims)], dims)
end
zeros(::Type{Obs}, dims::Tuple{}) where {(Obs<:MCObservable)} = reshape(Obs[Obs()], dims)
zeros(::Type{Obs}, dims::Integer...) where {(Obs<:MCObservable)} = zeros(Obs, dims)

function show(io::IO, obs::MCObservable)
    if !isempty(obs)
        print(io, mean(obs), " +/- ", stderror(obs))
    else
        print(io, "No Entries")
    end
end

function dump_plot(io::IO, obs::MCObservable; put_following_space::Bool=false,
                   newline::Bool=false)
    print(io, mean(obs), " ", stderror(obs))
    if put_following_space
        print(io, " ")
    end
    if newline
        println(io)
    end
end

const confidence_rate_1sigma = 0.5erf(0.5sqrt(2.0))

include("util.jl")
include("observableset.jl")
include("parsesigma.jl")
include("tiny.jl")
include("tinyvector.jl")
include("simple.jl")
include("simplevector.jl")
include("binning.jl")
include("binningvector.jl")
include("jackknife.jl")
include("jackknifevector.jl")

## The value type is constrained; a `<<(::MCObservable, ::Any)` catch-all is
## ambiguous with the `<<` methods that Base defines for its own types.
<<(obs::ScalarObservable, x::Real) = push!(obs, x)
<<(obs::VectorObservable, x::AbstractVector) = push!(obs, x)
