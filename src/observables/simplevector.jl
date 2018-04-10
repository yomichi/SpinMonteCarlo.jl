export SimpleVectorObservable, stddev

mutable struct SimpleVectorObservable <: VectorObservable
    bins :: Vector{Vector{Float64}}
    num :: Int64
    sum :: Vector{Float64}
    sum2 :: Vector{Float64}
end

SimpleVectorObservable() = SimpleVectorObservable(Vector{Float64}[], 0, Float64[], Float64[])

function reset!(obs :: SimpleVectorObservable)
    obs.bins = Vector{Vector{Float64}}[]
    obs.num = 0
    obs.sum = Float64[]
    obs.sum2 = Float64[]
    return obs
end

count(obs::SimpleVectorObservable) = obs.num

function push!(obs :: SimpleVectorObservable, value::Vector)
    if obs.num == 0
        obs.num = 1
        push!(obs.bins, value)
        obs.sum = deepcopy(value)
        obs.sum2 = squared(value)
    else
        push!(obs.bins, value)
        obs.num += 1
        for i in 1:length(obs.bins[1])
            obs.sum[i] += value[i]
            obs.sum2[i] += squared(value[i])
        end
    end
    return obs
end

function mean(obs::SimpleVectorObservable)
    if obs.num > 0
        return obs.sum .* (1.0/obs.num)
    else
        return [NaN]
    end
end

function var(obs::SimpleVectorObservable)
    if obs.num  > 1
        v = (obs.sum2 .- (obs.sum.^2)./obs.num)./(obs.num-1)
        return map!(maxzero, v)
    else
        return fill(NaN, length(obs.sum))
    end
end
stddev(obs::SimpleVectorObservable) = sqrt.(var(obs))
stderror(obs::SimpleVectorObservable) = sqrt.(var(obs)./count(obs))
function confidence_interval(obs::SimpleVectorObservable, confidence_rate :: Real)
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

function confidence_interval(obs::SimpleVectorObservable, confidence_rate_symbol::Symbol = :sigma1)
    n = parsesigma(confidence_rate_symbol)
    return confidence_interval(obs, erf(0.5n*sqrt(2.0)))
end

function merge!(obs::SimpleVectorObservable, other::SimpleVectorObservable)
    append!(obs.bins, other.bins)
    obs.num += other.num
    for i in 1:length(obs.bins[i])
        obs.sum[i] += other.sum[i]
        obs.sum2[i] += other.sum2[i]
    end
    return obs
end
merge(lhs::SimpleVectorObservable, rhs::SimpleVectorObservable) = merge!(deepcopy(lhs), rhs)

export SimpleVectorObservableSet
const SimpleVectorObservableSet =  MCObservableSet{SimpleVectorObservable}

function merge!(obs::SimpleVectorObservableSet, other::SimpleVectorObservableSet)
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
merge(lhs::SimpleVectorObservableSet, rhs::SimpleVectorObservableSet) = merge!(deepcopy(lhs), rhs)
