export JackknifeVector, JackknifeVectorSet

type JackknifeVector <: VectorObservable
  xs :: Vector{Vector{Float64}}
end

JackknifeVector() = JackknifeVector(Vector{Float64}[])

function JackknifeVector(jk::JackknifeVector, f::Function) 
  xs = similar(jk.xs)
  for i in 1:length(jk.xs)
    xs[i] = f(jk.xs[i])
  end
  JackknifeVector(xs)
end

function JackknifeVector(o::VectorObservable)
  if isempty(o)
    return JackknifeVector()
  else
    return JackknifeVector(jk_helper(o.bins))
  end
end

count(jk::JackknifeVector) = length(jk.xs)

function mean(jk::JackknifeVector)
  if isempty(jk) 
    return NaN
  else
    return mean(jk.xs)
  end
end
function stderror(jk::JackknifeVector)
  n = count(jk)
  if n == 0
    return NaN
  elseif n == 1
    return Inf
  else
    m2 = sumabs2(jk.xs)
    m2 /= n
    m = mean(jk)
    sigma2 = m2 - m.*m
    sigma2 *= n-1
    map!(maxzero, sigma2)
    return sqrt(sigma2)
  end
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
  @eval ($op)(jk::JackknifeVector) = JackknifeVector(jk, $op)
end

binary_functions = (
  :+, :-, :*, :/, :\
)

for op in ( :+, :-, :.+, :.- )
  eval( Expr(:import, :Base, op) )
  eval( Expr(:export, op) )
  @eval ($op)(jk::JackknifeVector, rhs::Real) = JackknifeVector(jk, lhs->($op)(lhs,rhs))
  @eval ($op)(jk::JackknifeVector, rhs::Vector) = JackknifeVector(jk, lhs->($op)(lhs,rhs))
  @eval ($op)(lhs::Real, jk::JackknifeVector) = JackknifeVector(jk, rhs->($op)(lhs,rhs))
  @eval ($op)(lhs::Vector, jk::JackknifeVector) = JackknifeVector(jk, rhs->($op)(lhs,rhs))
  @eval ($op)(lhs::JackknifeVector, rhs::JackknifeVector) = JackknifeVector( ($op)(lhs.xs, rhs.xs))
end
for op in ( :*, :/, :\)
  eval( Expr(:import, :Base, op) )
  eval( Expr(:export, op) )
  @eval ($op)(lhs::Real, jk::JackknifeVector) = JackknifeVector(jk, rhs->($op)(lhs,rhs))
  @eval ($op)(jk::JackknifeVector, rhs::Real) = JackknifeVector(jk, lhs->($op)(lhs,rhs))
end
for op in ( :.*, :./, :.\, :.^)
  eval( Expr(:import, :Base, op) )
  eval( Expr(:export, op) )
  @eval ($op)(lhs::Real, jk::JackknifeVector) = JackknifeVector(jk, rhs->($op)(lhs,rhs))
  @eval ($op)(lhs::Vector, jk::JackknifeVector) = JackknifeVector(jk, rhs->($op)(lhs,rhs))
  @eval ($op)(jk::JackknifeVector, rhs::Real) = JackknifeVector(jk, lhs->($op)(lhs,rhs))
  @eval ($op)(jk::JackknifeVector, rhs::Vector) = JackknifeVector(jk, lhs->($op)(lhs,rhs))
  @eval ($op)(lhs::JackknifeVector, rhs::JackknifeVector) = begin
    xs = similar(lhs.xs)
    for i in 1:length(lhs.xs)
      xs[i] = ($op)(lhs.xs[i], rhs.xs[i])
    end
    return JackknifeVector(xs)
  end
end

typealias JackknifeSet MCObservableSet{Jackknife}

jackknife(obs::VectorObservable) = JackknifeVector(obs)
function jackknife{Obs<:VectorObservable}(obsset :: MCObservableSet{Obs})
  JK = JackknifeSet()
  for (k,v) in obsset
    JK[k] = JackknifeVector(v)
  end
  return JK
end

