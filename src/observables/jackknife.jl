export Jackknife, JackknifeSet
export jackknife

type Jackknife <: ScalarObservable
    xs :: Vector{Float64}
end

Jackknife(jk::Jackknife, f::Function) = Jackknife(map(f,jk.xs))

function jk_helper(xs)
    s = sum(xs)
    n = length(xs)-1
    ret = similar(xs)
    for i in 1:n+1
        ret[i] = (s-xs[i])/n
    end
    return ret
end

function Jackknife(o::ScalarObservable)
    return count(o) > 0 ? Jackknife(jk_helper(o.bins)) : Jackknife(zeros(0))
end

count(jk::Jackknife) = length(jk.xs)

mean(jk::Jackknife) = mean(jk.xs)
function stderror(jk::Jackknife)
    n = count(jk)
    if n < 2
        return NaN
    else
        m2 = sum(abs2, jk.xs)
        m2 /= n
        m = mean(jk)
        sigma2 = m2 - m*m
        sigma2 *= n-1
        sigma2 = maxzero(sigma2)
        return sqrt(sigma2)
    end
end
function var(jk::Jackknife)
    n = count(jk)
    if n < 2
        return NaN
    else
        m = mean(jk)
        return sum(abs2, jk.xs.-m)*(n-1)
    end
end
stddev(jk::Jackknife) = sqrt(var(jk))

function confidence_interval(jk::Jackknife, confidence_rate::Real)
    q = 0.5+0.5*confidence_rate
    correction = quantile( TDist(count(jk)-1), q)
    serr = stderror(jk)
    return correction * serr
end
function confidence_interval(jk::Jackknife, confidence_rate_symbol::Symbol = :sigma1)
    n = parsesigma(confidence_rate_symbol)
    return confidence_interval(jk, erf(0.5n*sqrt(2.0)))
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
                   :erf, :erfc, :erfcx,
                   :erfinv, :erfcinv,
                   :gamma, :lgamma, :lfact,
                   :digamma, :invdigamma, :trigamma,
                   :airyai, :airyprime, :airyaiprime,
                   :airybi, :airybiprime,
                   :besselj0, :besselj1, 
                   :bessely0, :bessely1,
                   :eta, :zeta
                  )

for op in unary_functions
    eval( Expr(:import, :Base, op) )
    eval( Expr(:export, op) )
    @eval ($op)(jk::Jackknife) = Jackknife(jk, $op)
end

binary_functions = (
                    :+, :-, :*, :/, :\
                   )

for op in binary_functions
    eval( Expr(:import, :Base, op) )
    eval( Expr(:export, op) )
    @eval ($op)(jk::Jackknife, rhs::Real) = Jackknife(jk, lhs->($op)(lhs,rhs))
    @eval ($op)(lhs::Real, jk::Jackknife) = Jackknife(jk, rhs->($op)(lhs,rhs))
    op_bw = Symbol("."*string(op))
    @eval ($op)(lhs::Jackknife, rhs::Jackknife) = Jackknife( ($op_bw)(lhs.xs, rhs.xs))
end

const JackknifeSet = MCObservableSet{Jackknife}

jackknife(obs::ScalarObservable) = Jackknife(obs)
function jackknife(obsset :: MCObservableSet)
    JK = JackknifeSet()
    for (k,v) in obsset
        JK[k] = Jackknife(v)
    end
    return JK
end

