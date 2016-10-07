squared(x::Real) = x*x
squared(x::Vector) = x.*x
maxzero(x::Real) = max(x,zero(x))
