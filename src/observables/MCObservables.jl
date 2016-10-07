import Base: show, <<, push!, mean, var, count, isempty, merge, merge!, zero, zeros, sum
export MCObservable, ScalarObservable, VectorObservable
export mean, var, stderror, confidence_interval
export show, dump_plot
export merge, merge!

export confidence_interval

using Distributions

abstract MCObservable
abstract ScalarObservable <: MCObservable
abstract VectorObservable <: MCObservable

isempty(obs::MCObservable) = count(obs) == 0

zero{Obs<:MCObservable}(::Type{Obs}) = Obs()
zero(obs::MCObservable) = zero(typeof(obs))
function zeros{Obs<:MCObservable}(::Type{Obs}, dim...)
  reshape(Obs[Obs() for i in 1:prod(dim)], dim)
end

function show(io::IO, obs::MCObservable)
  if !isempty(obs)
    print(io, mean(obs), " +/- ", stderror(obs))
  else
    print(io, "No Entries")
  end
end

function dump_plot(io::IO, obs::MCObservable)
  print(io, mean(obs), " ", stderror(obs))
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

## these three definitions are needed to resolve ambiguousness with Base.<<(Any,Integer)
<<(obs::MCObservable, x::Int32) = push!(obs,x)
<<(obs::MCObservable, x::Int64) = push!(obs,x)
<<(obs::MCObservable, x::Integer) = push!(obs,x)

<<(obs::MCObservable, x) = push!(obs,x)

