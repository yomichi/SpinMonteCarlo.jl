export DoS
export learnEE
export measureEE

mutable struct DoS{T<:Real}
    log_g :: Vector{Float64}
    hist :: Vector{UInt128}
    lower_bounds :: T
    upper_bounds :: T
    nbin :: UInt64
    delta :: T
end
function DoS(lower_bounds::T, upper_bounds::T, delta::T) where {T<:Real}
    nbin = (upper_bounds - lower_bounds) / delta
    return DoS{T}(zeros(nbin), zeros(UInt128), lower_bounds, upper_bounds, nbin, delta)
end
function DoS(lower_bounds::T, upper_bounds::T, delta::T) where {T<:Integer}
    nbin = (upper_bounds - lower_bounds) ÷ delta + 1
    return DoS{T}(zeros(nbin), zeros(UInt128, nbin), lower_bounds, upper_bounds, nbin, delta)
end

valuetype(dos::DoS{T}) where {T<:Real} = T

function value2index(dos::DoS, value::Real; check_range::Bool = true)
    if check_range
        if ! ( dos.lower_bounds <= value <= dos.upper_bounds )
            error()
        end
    else
        value = clamp(value, dos.lower_bounds, dos.upper_bounds)
    end
    value = ifelse(value == dos.upper_bounds, prevfloat(value), value)
    return ceil(UInt64, (value - dos.lower_bounds) / dos.delta)
end

function value2index(dos::DoS{T}, value::T; check_range::Bool = true) where {T<:Integer}
    if check_range
        if ! ( dos.lower_bounds <= value <= dos.upper_bounds )
            error()
        end
    else
        value = clamp(value, dos.lower_bounds, dos.upper_bounds)
    end
    return (value - dos.lower_bounds)÷dos.delta + 1
end

function index2value(dos::DoS, index::Integer)
    return dos.lower_bounds + (index-0.5) * dos.delta
end

function index2value(dos::DoS{T}, index::Integer) where {T<:Integer}
    return dos.lower_bounds + (index-1)*dos.delta
end

function visit!(dos::DoS, value::Real, α::Real=0.0; check_range::Bool = true)
    index = value2index(dos, value, check_range=check_range)
    dos.log_g[index] += α
    dos.hist[index] += 1
end

function doshist(dos::DoS, value::Real; check_range::Bool = true)
    index = value2index(dos, value, check_range=check_range)
    return dos.log_g[index], dos.hist[index]
end

@inline log_g(dos::DoS, value::Real; check_range::Bool = true) = doshist(dos, value, check_range=check_range)[1]
@inline hist(dos::DoS, value::Real; check_range::Bool = true) = doshist(dos, value, check_range=check_range)[2]

function check_flat(dos::DoS, α; criteria::Real = 0.8, nonzero_ratio::Real = 0.5, verbose::Bool=true)
    hmin = typemax(UInt128)
    hmax = typemin(UInt128)
    imin = 0
    imax = 0
    hsum = zero(UInt128)
    hnum = 0
    for (i, h) in enumerate(dos.hist)
        if h > hmax
            hmax = h
            imax = i
        end
        if !iszero(h)
            if h < hmin
                hmin = h
                imin = i
            end
            hsum += h
            hnum += 1
        end
    end
    if hnum == 0
        return false
    end
    hmean = hsum / hnum
    flatness = hmin / hmax
    if hnum < nonzero_ratio * dos.nbin
        if verbose
            @printf("Too small region is explored (%.3f < %.3f) [hmax = %d@%d, hmin = %d@%d, α = %g] \n",
                    hnum/dos.nbin, nonzero_ratio, hmax, imax, hmin, imin, α)
        end
        return false
    else
        if flatness < criteria
            if verbose
                @printf("Not flat (%.3f < %.3f)  [hmax = %d@%d, hmin = %d@%d, α = %g] \n",
                        flatness, criteria, hmax, imax, hmin, imin, α)
            end
            return false
        else
            @printf("Flat     (%.3f >= %.3f) [hmax = %d@%d, hmin = %d@%d, α = %g] \n",
                    flatness, criteria, hmax, imax, hmin, imin, α)
            return true
        end
    end
end

function reset_hist!(dos::DoS)
    dos.hist .= zero(UInt128)
    return dos
end

function renormalize!(dos::DoS)
    dos.log_g .-= median(dos.log_g)
    return dos
end

function updateEE!(model::Model, present, dos::DoS, dosobsname, p)
    action = nextaction(model)
    diff = localchange(model, action, p...)[dosobsname]
    if rand() < exp(log_g(dos, present) - log_g(dos, present+diff))
        present += diff
        accept!(model, action)
    end
    return present
end


function learnEE(param::Parameter)
    model = param["Model"](param)
    return learnEE(model, param)
end

function learnEE(model::Model, param::Parameter)
    dos = init_extended_ensemble(model, param)
    dosobsname = param["Observable for Extended Ensemble"]

    p = convert_parameter(model, param)

    α = get(param, "WL Start Refinement Factor", 1.0)
    α_last = get(param, "WL Minimum Refinement Factor", 1e-4)

    criteria = get(param, "WL Criteria for Flat", 0.8)
    nonzero_ratio = get(param, "WL Nonzero Ratio", 0.4)

    interval_calc_full_state = 100

    istage = 1

    while α >= α_last
        reset_hist!(dos)
        iter = 0
        present = simple_estimator(model, p..., nothing)[dosobsname]
        while true
            for i in 1:numsites(model)
                present = updateEE!(model, present, dos, dosobsname, p)
                visit!(dos, present, α)
            end
            if check_flat(dos, α, criteria = criteria, nonzero_ratio = nonzero_ratio, verbose=true)
                break
            end
            iter += 1
            if iter % interval_calc_full_state == 0
                present = simple_estimator(model, p..., nothing)[dosobsname]
                iter = 0
            end
        end
        renormalize!(dos)
        open("dos_$(istage).dat", "w") do io
            if valuetype(dos) <: Integer
                for (i, (g, h)) in enumerate(zip(dos.log_g, dos.hist))
                    v = index2value(dos, i)
                    @printf(io, "%d %.15f %d\n", v, g, h)
                end
            else
                for (i, (g, h)) in enumerate(zip(dos.log_g, dos.hist))
                    @printf(io, "%.15f %.15f %d\n", v, g, h)
                end
            end
        end
        α *= 0.5
        istage += 1
    end
end

function measureEE(model::Model, dos::DoS, param::Parameter)
    MCS = get(param, "MCS", 65536)
    Thermalization = get(param, "Thermalization", MCS >> 3)

    dosobsname = param["Observable for Extended Ensemble"]
    p = convert_parameter(model, param)
    present = simple_estimator(model, p..., nothing)[dosobsname]

    for mcs in 1:Thermalization
        for i in 1:numsites(model)
            present = updateEE!(model, present, dos, dosobsname, p)
        end
        present = simple_estimator(model, p..., nothing)[dosobsname]
    end

    obs = BinningObservableSet()

    for mcs in 1:MCS
        for i in 1:numsites(model)
            present = updateEE!(model, present, dos, dosobsname, p)
        end
        lobs = simple_estimator(model, p..., nothing)
        lobs["Identity"] = 1.0
        present = lobs[dosobsname]
        accumulateObservablesEE!(model, log_g(dos, present), obs, lobs)
    end
    
    jk = jackknife(obs)
    for key in keys(jk)
        if ! startswith(key, "Wg ")
            continue
        end
        obsname = string(key[4:end])
        jk[obsname] = jk[key] / jk["Wg Identity"]
    end

    jk = postproc(model, param, jk)

    return jk
end

@doc """
    accumulateObservablesEE!(model, logg, obs::MCObservableSet, localobs::Dict)

Accumulates `localobs` into `obs`. For example, `obs["Energy"] << localobs["Energy"]`.
"""
function accumulateObservablesEE!(::Model, logg, obs::MCObservableSet, localobs::Measurement)
    if length(obs) < 3
        @inbounds for key in keys(localobs)
            obsname = "Wg " * key
            makeMCObservable!(obs, obsname)
            v = localobs[key]
            v *= exp(localobs["Log Boltzmann Weight"] + logg)
            obs[obsname] << v
        end
    else
        @inbounds for key in keys(localobs)
            v = localobs[key]
            v *= exp(localobs["Log Boltzmann Weight"] + logg)
            obs["Wg " * key] << v
        end
    end
    return obs
end
