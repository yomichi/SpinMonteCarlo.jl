@inline squared(x::Real) = x*x
@inline squared(x::Vector) = x.*x
@inline maxzero(x::Real) = max(x,zero(x))
