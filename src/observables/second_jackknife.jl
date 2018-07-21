export SecondJackknife, SecondJackknifeSet
export second_jackknife

mutable struct SecondJackknife <: MCObservable
    x_i  :: Vector{Float64}
    x_ij :: Matrix{Float64}
    b_i  :: Vector{Float64}
end

function SecondJackknife(x_ij :: Matrix{Float64})
    n = size(x_ij,1)
    x_i = zeros[n]
    b_i = zeros[n]
    for i in 1:n
        x_i[i] = x_ij[i,i]
        b_i[i] = x_ij[i,i]
        b_i[i] -= sum(x_ij[1:i-1,i])
        b_i[i] -= sum(x_ij[i,i+1:n])
    end
    return SecondJackknife(x_i, x_ij, b_i)
end

SecondJackknife(jk::SecondJackknife, f::Function) = SecondJackknife(map(f,jk.x_ij))

function SecondJackknife(b::BinningObservable)
    if count(b) > 0
        level = maxlevel(b)
        xs = b.bins[level]
        if b.nincomplete[level] > 0
            push!(xs, b.incompletes[level]/b.nincomplete[level])
        end
        s = sum(xs)
        n = length(xs)
        n1inv = 1.0/(n-1)
        n2inv = 1.0/(n-2)
        jk_ij = zeros(n,n)
        for i in 1:n
            xi = xs[i]
            jk_ij[i,i] = (s-xi)*n1inv
            for j in i+1:n
                jk_ij[i,j] = (s-xi-xs[j])*n2inv
            end
        end
        return SecondJackknife(jk_ij)
    else
        return Jackknife(zeros(0,0))
    end
end

count(jk::SecondJackknife) = length(jk.x_ij)
isempty(jk::SecondJackknife) = isempty(jk.x_ij)

function mean(jk::SecondJackknife)
    if isempty(jk) 
        throw(DomainError())
    else
        return mean(jk.x_i + jk.b_i)
    end
end
function stderr(jk::SecondJackknife)
    n = count(jk)
    if n == 0
        throw(DomainError())
    elseif n == 1
        return Inf
    else
        sums = mapreduce(x->[x,x*x], +, jk.x_i + jk.b_i)
        sums /= n
        sigma2 = sums[2] - sums[1]*sums[1]
        sigma2 *= n-1
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
                  )

for op in unary_functions
    eval( Expr(:import, :Base, op) )
    eval( Expr(:export, op) )
    @eval ($op)(jk::SecondJackknife) = SecondJackknife(jk, $op)
end

binary_functions = (
                    :+, :-, :*, :/, :\
                   )

for op in binary_functions
    eval( Expr(:import, :Base, op) )
    eval( Expr(:export, op) )
    @eval ($op)(jk::SecondJackknife, rhs::Real) = SecondJackknife(jk, lhs->($op)(lhs,rhs))
    @eval ($op)(lhs::Real, jk::SecondJackknife) = SecondJackknife(jk, rhs->($op)(lhs,rhs))
    op_bw = symbol("."*string(op))
    @eval ($op)(lhs::SecondJackknife, rhs::SecondJackknife) = SecondJackknife( ($op_bw)(lhs.xs, rhs.xs))
end

const SecondJackknifeSet = MCObservableSet{SecondJackknife}

function second_jackknife(obsset :: BinningObservableSet)
    JK = SecondJackknifeSet()
    for (k,v) in obsset
        JK[k] = SecondJackknife(v)
    end
    return JK
end

