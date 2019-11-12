export DoS
export visit!, logg, hist, logg_hist, reset_hist!, renormalize!
export copy_bins!

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

function value2ratio(dos::DoS, value)
    return (value - dos.lower_bounds) / (dos.upper_bounds - dos.lower_bounds)
end

index2ratio(dos::DoS, index) = value2ratio(dos, index2value(dos, index))

function visit!(dos::DoS, value::Real, α::Real=0.0; check_range::Bool = true)
    index = value2index(dos, value, check_range=check_range)
    dos.log_g[index] += α
    dos.hist[index] += 1
end

function logg_hist(dos::DoS, value::Real; check_range::Bool = true)
    index = value2index(dos, value, check_range=check_range)
    return dos.log_g[index], dos.hist[index]
end

@inline logg(dos::DoS, value::Real; check_range::Bool = true) = logg_hist(dos, value, check_range=check_range)[1]
@inline hist(dos::DoS, value::Real; check_range::Bool = true) = logg_hist(dos, value, check_range=check_range)[2]

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

function dump_dos(filename, dos::DoS; skip_empty::Bool=true)
    open(filename, "w") do io
        if SpinMonteCarlo.valuetype(dos) <: Integer
            for (i, (g, h)) in enumerate(zip(dos.log_g, dos.hist))
                if skip_empty && h > 0
                    v = SpinMonteCarlo.index2value(dos, i)
                    @printf(io, "%d %.15f %d\n", v, g, h)
                end
            end
        else
            for (i, (g, h)) in enumerate(zip(dos.log_g, dos.hist))
                if skip_empty && h > 0
                    v = SpinMonteCarlo.index2value(dos, i)
                    @printf(io, "%.15f %.15f %d\n", v, g, h)
                end
            end
        end
    end
end


function save_dos(filename, dos::DoS)
    jldopen(filename, "w") do file
        file["log g"] = dos.log_g
        file["hist"] = dos.hist
        file["lower bounds"] = dos.lower_bounds
        file["upper bounds"] = dos.upper_bounds
        file["nbin"] = dos.nbin
        file["delta"] = dos.delta
    end
end

function load_dos(filename)
    return jldopen(filename) do file
        log_g = file["log g"]
        hist = file["hist"]
        lb = file["lower bounds"]
        ub = file["upper bounds"]
        nbin = file["nbin"]
        d = file["delta"]
        dos = DoS(promote(lb, ub, d)...)
        dos.log_g .= log_g
        dos.hist .= hist
        return dos
    end
end

function copy_dos!(dst::DoS, src::DoS)
    for idst in 1:dst.nbin
        r = (idst-1) / (dst.nbin-1)
        isrc_re = r * (src.nbin-1) + 1.0
        if floor(isrc_re) == isrc_re
            isrc_re = ifelse(isrc_re > 1.0, prevfloat(isrc_re), nextfloat(isrc_re))
        end
        isrc_floor = floor(Int, isrc_re)
        isrc = isrc_floor
        while isrc_floor >= 1 && src.hist[isrc_floor] == 0
            isrc_floor -= 1
        end
        isrc_ceil = ceil(Int, isrc_re)
        while isrc_ceil <= src.nbin && src.hist[isrc_ceil] == 0
            isrc_ceil += 1
        end

        dst.log_g[idst] = src.log_g[isrc_floor] + (isrc_re - isrc_floor) * (src.log_g[isrc_ceil] - src.log_g[isrc_floor]) / (isrc_ceil - isrc_floor)
        dst.hist[idst] = floor(UInt128, src.hist[isrc_floor] + (isrc_re - isrc_floor) * (src.hist[isrc_ceil] - src.hist[isrc_floor]) / (isrc_ceil - isrc_floor))
    end
    return dst
end
