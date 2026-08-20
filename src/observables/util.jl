@inline squared(x::Real) = x * x
@inline squared(x::Vector) = x .* x
@inline maxzero(x::Real) = max(x, zero(x))

"""
    linear_intercept(xs, ys)

Fit `y = a + b*x` by unweighted ordinary least squares and return
`(a, stderror_of_a)`. The error is a 1-sigma standard error.

Throws `ArgumentError` if fewer than three points are given, since the
standard error needs at least one degree of freedom (`dof = n - 2`), or if
any input is non-finite.
"""
function linear_intercept(xs::AbstractVector{<:Real}, ys::AbstractVector{<:Real})
    n = length(xs)
    n == length(ys) ||
        throw(ArgumentError("xs and ys must have the same length, got $n and $(length(ys))"))
    n >= 3 ||
        throw(ArgumentError("linear extrapolation needs at least 3 points, got $n"))

    for (i, y) in enumerate(ys)
        isfinite(y) ||
            throw(ArgumentError("non-finite value $y at index $i; a fit cannot be performed (this usually means a bin level holds fewer than two bins)"))
    end
    for (i, x) in enumerate(xs)
        isfinite(x) || throw(ArgumentError("non-finite x value $x at index $i"))
    end

    xbar = sum(xs) / n
    ybar = sum(ys) / n
    Sxx = sum(squared(x - xbar) for x in xs)
    Sxx > 0 || throw(ArgumentError("all x values are identical, cannot fit a line"))
    Sxy = sum((x - xbar) * (y - ybar) for (x, y) in zip(xs, ys))

    slope = Sxy / Sxx
    intercept = ybar - slope * xbar

    dof = n - 2
    rss = sum(squared(y - (intercept + slope * x)) for (x, y) in zip(xs, ys))
    varhat = rss / dof
    stderr_intercept = sqrt(varhat * (1 / n + squared(xbar) / Sxx))

    return intercept, stderr_intercept
end
