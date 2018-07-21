export TinyVectorObservable, stddev

mutable struct TinyVectorObservable <: VectorObservable
    num :: Int64
    sum :: Vector{Float64}
    sum2 :: Vector{Float64}
end

TinyVectorObservable() = TinyVectorObservable(0, Float64[], Float64[])
zero(::Type{TinyVectorObservable}) = TinyVectorObservable()
zero(o::TinyVectorObservable) = TinyVectorObservable()
zeros(::Type{TinyVectorObservable}, dims...) = reshape([zero(TinyVectorObservable) for i in 1:prod(dims)],dims)

function reset!(obs :: TinyVectorObservable)
    obs.num = 0
    obs.sum = Float64[]
    obs.sum2 = Float64[]
    return obs
end

count(obs::TinyVectorObservable) = obs.num

function push!(obs :: TinyVectorObservable, value::Vector)
    if obs.num == 0
        obs.num = 1
        obs.sum = deepcopy(value)
        obs.sum2 = squared(value)
    else
        obs.num += 1
        for i in 1:length(obs.sum)
            obs.sum[i] += value[i]
            obs.sum2[i] += squared(value[i])
        end
    end
    return obs
end

function mean(obs::TinyVectorObservable)
    if obs.num > 0
        return obs.sum .* (1.0/obs.num)
    else
        return [NaN]
    end
end

function var(obs::TinyVectorObservable)
    if obs.num  > 1
        v = (obs.sum2 .- (obs.sum.^2)./obs.num)./(obs.num-1)
        return map!(maxzero, v)
    else
        return fill(NaN, length(obs.sum))
    end
end
stddev(obs::TinyVectorObservable) = sqrt.(var(obs))
stderror(obs::TinyVectorObservable) = sqrt.(var(obs)./count(obs))
function confidence_interval(obs::TinyVectorObservable, confidence_rate :: Real)
    if count(obs) == 0
        return [Inf]
    elseif count(obs) == 1
        return fill(Inf, length(obs.sum))
    end
    q = 0.5 + 0.5confidence_rate
    correction = quantile( TDist(obs.num - 1), q)
    serr = stderror(obs)
    return correction * serr
end

function confidence_interval(obs::TinyVectorObservable, confidence_rate_symbol::Symbol = :sigma1)
    n = parsesigma(confidence_rate_symbol)
    return confidence_interval(obs, erf(0.5n*sqrt(2.0)))
end

function merge!(obs::TinyVectorObservable, other::TinyVectorObservable)
    obs.num += other.num
    for i in 1:length(obs.sum)
        obs.sum[i] += other.sum[i]
        obs.sum2[i] += other.sum2[i]
    end
    return obs
end
merge(lhs::TinyVectorObservable, rhs::TinyVectorObservable) = merge!(deepcopy(lhs), rhs)

export TinyVectorObservableSet
const TinyVectorObservableSet = MCObservableSet{TinyVectorObservable}

function merge!(obs::TinyVectorObservableSet, other::TinyVectorObservableSet)
    obs_names = Set(keys(obs))
    union!(obs_names, Set(keys(other)))
    for name in obs_names
        if !haskey(obs, name)
            # in other only
            obs[name] = deepcopy(other[name])
        elseif haskey(other, name)
            # in both
            merge!(obs[name], other[name])

            # else
            # in obs only
            # NOTHING to do
        end
    end
    return obs
end
merge(lhs::TinyVectorObservableSet, rhs::TinyVectorObservableSet) = merge!(deepcopy(lhs), rhs)
