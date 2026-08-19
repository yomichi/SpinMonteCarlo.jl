function parsesigma(confidence_rate_symbol::Symbol=:sigma1)
    m = match(r"^sigma(\d+)$", string(confidence_rate_symbol))
    n = m === nothing ? 0 : something(tryparse(Int, m.captures[1]), 0)
    if n <= 0
        throw(ErrorException("invalid parameter confidence_rate should be :sigmaN, where N is a positive integer"))
    end
    return n
end
