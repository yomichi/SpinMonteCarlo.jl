export JackknifeVector, JackknifeVectorSet

mutable struct JackknifeVector <: VectorObservable
    xs :: Vector{Vector{Float64}}
end

function jk_helper(xs::Vector{Vector{Float64}})
    ndata = length(xs)
    nobs = length(xs[1])
    ret = Vector{Float64}[zeros(nobs) for i in 1:ndata]

    s = sum(xs)
    coeff = 1.0 / (ndata - 1)
    n = length(xs)-1
    ret = similar(xs)
    for i in 1:n+1
        ret[i] = coeff .* (s.-xs[i])
    end
    return ret
end

JackknifeVector() = JackknifeVector(Vector{Float64}[])

function JackknifeVector(f::Function, jks::JackknifeVector...)
    if isempty(jks)
        return JackknifeVector()
    end
    xs = similar(jks[1].xs)
    for i in 1:length(jks[1].xs)
        arg = [jk.xs[i] for jk in jks]
        xs[i] = broadcast(f, arg...)
    end
    JackknifeVector(xs)
end

function JackknifeVector(jk::JackknifeVector, f::Function) 
    Base.depwarn("JackknifeVector(jk::JackknifeVector, f::Function) is deprecated. Use JackknifeVector(f, jk) instead.", :JackknifeVector)
    return JackknifeVector(f, jk)
end

function JackknifeVector(f::Function, lhs::JackknifeVector, rhs::Jackknife)
    nbins = count(lhs)
    nobs = length(lhs.xs[1])
    xs = [zeros(nobs) for _ in lhs.xs]
    for i in 1:nbins
        arg = [lhs.xs[i], rhs.xs[i]]
        xs[i] .= broadcast(f, arg...)
    end
    return JackknifeVector(xs)
end

function JackknifeVector(f::Function, lhs::Jackknife, rhs::JackknifeVector)
    nbins = count(lhs)
    nobs = length(rhs.xs[1])
    xs = [zeros(nobs) for _ in rhs.xs]
    for i in 1:nbins
        arg = [lhs.xs[i], rhs.xs[i]]
        xs[i] .= broadcast(f, arg...)
    end
    return JackknifeVector(xs)
end


function JackknifeVector(o::VectorObservable)
    if isempty(o)
        return JackknifeVector()
    else
        return JackknifeVector(jk_helper(o.bins))
    end
end
function JackknifeVector(b::BinningVectorObservable)
    if isempty(b)
        return JackknifeVector()
    else
        if b.lastbin == b.binsize
            return JackknifeVector(jk_helper(b.bins))
        else
            bins = [bs[1:end-1] for bs in b.bins]
            return JackknifeVector(jk_helper(bins))
        end
    end
end

count(jk::JackknifeVector) = length(jk.xs)

function mean(jk::JackknifeVector)
    if isempty(jk) 
        return [NaN]
    else
        res = zeros(length(jk.xs[1]))
        for data in jk.xs
            res .+= data
        end
        return res .* (1.0 / length(jk.xs))
    end
end
function var(jk::JackknifeVector)
    n = count(jk)
    if n < 2
        return [NaN]
    else
        m = mean(jk)
        s = similar(m)
        for data in jk.xs
            s .+= (data .- m) .^ 2
        end
        s .*= (n-1.0)/n
        return s
    end
end
stddev(jk::JackknifeVector) = sqrt.(var(jk))
function stderror(jk::JackknifeVector)
    # n = count(jk)
    # if n == 0
    #     return NaN
    # elseif n == 1
    #     return Inf
    # else
    #     m2 = sum(abs2, jk.xs)
    #     m2 /= n
    #     m = mean(jk)
    #     sigma2 = m2 - m.*m
    #     sigma2 *= n-1
    #     map!(maxzero, sigma2)
    #     return sqrt(sigma2)
    # end
    return stddev(jk)
end


unary_functions = (
                   :-,
                   :sin, :cos, :tan,
                   :sind, :cosd, :tand,
                   :sinpi, :cospi,
                   :sinh, :cosh, :tanh,
                   :asin, :acos, :atan,
                   :asind, :acosd, :atand,
                   :sec, :csc, :cot,
                   :secd, :cscd, :cotd,
                   :asec, :acsc, :acot,
                   :asecd, :acscd, :acotd,
                   :sech, :csch, :coth,
                   :asinh, :acosh, :atanh,
                   :asech, :acsch, :acoth,
                   :sinc, :cosc,
                   :log, :log2, :log10, :log1p,
                   :exp, :exp2, :exp10, :expm1,
                   :abs, :abs2,
                   :sqrt, :cbrt,
                  )

for op in unary_functions
    @eval Base.$op(jk::JackknifeVector) = JackknifeVector($op, jk)
end

binary_functions = (
                    :+, :-, :*, :/, :\
                   )

import Base.broadcast
for op in ( :+, :- )
    @eval Base.$op(jk::JackknifeVector, rhs::Real) = JackknifeVector(lhs->($op)(lhs,rhs), jk)
    @eval Base.$op(jk::JackknifeVector, rhs::Vector) = JackknifeVector(lhs->($op)(lhs,rhs), jk)
    @eval Base.$op(lhs::Real, jk::JackknifeVector) = JackknifeVector(rhs->($op)(lhs,rhs), jk)
    @eval Base.$op(lhs::Vector, jk::JackknifeVector) = JackknifeVector(rhs->($op)(lhs,rhs), jk)
    @eval Base.$op(lhs::Jackknife, rhs::JackknifeVector) = JackknifeVector($op, lhs, rhs)
    @eval Base.$op(lhs::JackknifeVector, rhs::Jackknife) = JackknifeVector($op, lhs, rhs)
    @eval Base.$op(lhs::JackknifeVector, rhs::JackknifeVector) = JackknifeVector($op, lhs, rhs)
    @eval broadcast(::typeof($op), jk::JackknifeVector, rhs::Real) = JackknifeVector(lhs->($op)(lhs,rhs), jk)
    @eval broadcast(::typeof($op), jk::JackknifeVector, rhs::Vector) = JackknifeVector(lhs->($op)(lhs,rhs), jk)
    @eval broadcast(::typeof($op), lhs::Real, jk::JackknifeVector) = JackknifeVector(rhs->($op)(lhs,rhs), jk)
    @eval broadcast(::typeof($op), lhs::Vector, jk::JackknifeVector) = JackknifeVector(rhs->($op)(lhs,rhs), jk)
    @eval broadcast(::typeof($op), lhs::Jackknife, rhs::JackknifeVector) = JackknifeVector($op, lhs, rhs)
    @eval broadcast(::typeof($op), lhs::JackknifeVector, rhs::Jackknife) = JackknifeVector($op, lhs, rhs)
    @eval broadcast(::typeof($op), lhs::JackknifeVector, rhs::JackknifeVector) = JackknifeVector($op, lhs, rhs)
end
for op in ( :*, :/, :\)
    @eval Base.$op(lhs::Real, jk::JackknifeVector) = JackknifeVector(rhs->($op)(lhs,rhs), jk)
    @eval Base.$op(jk::JackknifeVector, rhs::Real) = JackknifeVector(lhs->($op)(lhs,rhs), jk)
    @eval Base.$op(lhs::JackknifeVector, rhs::Jackknife) = JackknifeVector($op, lhs, rhs)
    @eval Base.$op(lhs::Jackknife, rhs::JackknifeVector) = JackknifeVector($op, lhs, rhs)
end
for op in ( :*, :/, :\)
    @eval broadcast(::typeof($op), lhs::Real, jk::JackknifeVector) = JackknifeVector(rhs->($op)(lhs,rhs), jk)
    @eval broadcast(::typeof($op), lhs::Vector, jk::JackknifeVector) = JackknifeVector(rhs->($op)(lhs,rhs), jk)
    @eval broadcast(::typeof($op), jk::JackknifeVector, rhs::Real) = JackknifeVector(lhs->($op)(lhs,rhs))
    @eval broadcast(::typeof($op), jk::JackknifeVector, rhs::Vector) = JackknifeVector(jk, lhs->($op)(lhs,rhs))
    @eval broadcast(::typeof($op), lhs::Jackknife, rhs::JackknifeVector) = JackknifeVector($op, lhs, rhs)
    @eval broadcast(::typeof($op), lhs::JackknifeVector, rhs::Jackknife) = JackknifeVector($op, lhs, rhs)
    @eval broadcast(::typeof($op), lhs::JackknifeVector, rhs::JackknifeVector) = JackknifeVector($op, lhs, rhs)
end

const JackknifeVectorSet = MCObservableSet{JackknifeVector}

jackknife(obs::VectorObservable) = JackknifeVector(obs)
function jackknife(obsset :: MCObservableSet{Obs}) where (Obs<: VectorObservable)
    JK = JackknifeSet()
    for (k,v) in obsset
        JK[k] = jackknife(v)
    end
    return JK
end
