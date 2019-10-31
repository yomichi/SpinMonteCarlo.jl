mutable struct DoS
    log_g :: Vector{Float64}
    hist :: Vector{UInt128}
    lower_bounds :: Float64
    upper_bounds :: Float64
    nbin :: UInt64
    delta :: Float64
end
function DoS(lower_bounds, upper_bounds, nbin)
    return DoS(zeros(nbin), zeros(UInt128), lower_bounds, upper_bounds, nbin, (upper_bounds-lower_bounds)/(nbin-1))
end

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

function visit!(dos::DoS, value::Real, α::Real=0.0; check_range::Bool = true)
    index = value2index(dos, value, check_range=check_range)
    dos.dos[index] += α
    dos.hist[index] += 1
end

function doshist(dos::DoS, value::Real; check_range::Bool = true)
    index = value2index(dos, value, check_range=check_range)
    return dos.dos[index], dos.hist[index]
end

@inline log_g(dos::DoS, value::Real; check_range::Bool = true) = doshist(dos, value, check_range=check_range)[1]
@inline hist(dos::DoS, value::Real; check_range::Bool = true) = doshist(dos, value, check_range=check_range)[2]

function check_flat(dos::DoS; criteria::Real = 0.8, nonzero_ratio::Real = 0.5, verbose::Bool=true)
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
            @printf("Too small region is explored (%.3f < %.3f) \n", hnum/dos.nbin, nonzero_ratio)
        end
        return false
    else
        if flatness < criteria
            if verbose
                @printf("Not flat (%.3f < %.3f)  [hmax = %d@%d, hmin = %d@%d] \n", flatness, criteria, hmax, imax, hmin, imin)
            end
            return false
        else
            if verbose
                @printf("Flat     (%.3f >= %.3f) [hmax = %d@%d, hmin = %d@%d] \n", flatness, criteria, hmax, imax, hmin, imin)
            end
            return true
        end
    end
end

function reset_hist!(dos::DoS)
    dos.hist .= zero(UInt128)
end

function renormalize_logg!(dos::DoS)
    dos.log_g .-= median(dos.log_g)
end

function learnDoS(model::Model, param::Parameter)
    dosobsname = param["DoS Name"]

    p = convert_parameter(model, param)
    dos = DoS()

    present = simple_estimator(model, p..., nothing)[dosobsname]

    α = 1.0
    α_last = 1e-4

    criteria = 0.8
    nonzero_ratio = 0.4

    while α > α_last
        while !flatness(dos, criteria = criteria, nonzero_ratio = nonzero_ratio, verbose=true)
            action = nextaction(model)
            diff = localchange(model, action)
            if rand() < exp(log_g(dos, present) - log_g(dos, present+diff))
                present += diff
                accept!(model, action)
            end
            visit!(dos, present, α)
        end
        reset_hist!(dos)
        renormalize_logg!(dos)
        α *= 0.5
    end
end
