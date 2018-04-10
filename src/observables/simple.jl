export SimpleObservable, stddev

mutable struct SimpleObservable <: ScalarObservable
    bins :: Vector{Float64}
    num :: Int64
    sum :: Float64
    sum2 :: Float64
end

SimpleObservable() = SimpleObservable(zeros(0), 0, 0.0, 0.0)

function reset!(obs :: SimpleObservable)
    obs.bins = zeros(0)
    obs.num = 0
    obs.sum = 0.0
    obs.sum2 = 0.0
    return obs
end

count(obs::SimpleObservable) = obs.num

function push!(obs :: SimpleObservable, value) 
    push!(obs.bins, value)
    obs.num += 1
    obs.sum += value
    obs.sum2 += value^2
    return obs
end

function mean(obs::SimpleObservable)
    if obs.num > 0
        return obs.sum / obs.num
    else
        return NaN
    end
end

function var(obs::SimpleObservable)
    if obs.num  > 1
        v = (obs.sum2 - obs.sum*obs.sum/obs.num)/(obs.num-1)
        return maxzero(v)
    elseif obs.num < 2
        return NaN
    end
end
stddev(obs::SimpleObservable) = sqrt(var(obs))
stderror(obs::SimpleObservable) = sqrt(var(obs)/count(obs))
function confidence_interval(obs::SimpleObservable, confidence_rate :: Real)
    if count(obs) < 2
        return Inf
    end
    q = 0.5 + 0.5confidence_rate
    correction = quantile( TDist(obs.num - 1), q)
    serr = stderror(obs)
    return correction * serr
end

function confidence_interval(obs::SimpleObservable, confidence_rate_symbol::Symbol = :sigma1)
    n = parsesigma(confidence_rate_symbol)
    return confidence_interval(obs, erf(0.5n*sqrt(2.0)))
end

function merge!(obs::SimpleObservable, other::SimpleObservable)
    append!(obs.bins, other.bins)
    obs.num += other.num
    obs.sum += other.sum
    obs.sum2 += other.sum2
    return obs
end
merge(lhs::SimpleObservable, rhs::SimpleObservable) = merge!(deepcopy(lhs), rhs)

export SimpleObservableSet
const SimpleObservableSet = MCObservableSet{SimpleObservable}

function merge!(obs::SimpleObservableSet, other::SimpleObservableSet)
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
merge(lhs::SimpleObservableSet, rhs::SimpleObservableSet) = merge!(deepcopy(lhs), rhs)
