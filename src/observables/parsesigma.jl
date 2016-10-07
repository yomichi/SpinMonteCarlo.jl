function parsesigma(confidence_rate_symbol::Symbol = :sigma1)
  if confidence_rate_symbol == :sigma1
    n = 1
  elseif confidence_rate_symbol == :sigma2
    n = 2
  elseif confidence_rate_symbol == :sigma3
    n = 3
  else
    m = match( r"^sigma(\d+)$", string(confidence_rate_symbol) )
    if m == nothing
      n = 0
    end
    dstr = m.captured[1]
    n = try
      parseint(dstr)
    catch ex
      0
    end
  end
  if n <= 0
    throw(ErrorException("invalid parameter confidence_rate should be :sigmaN, where N is a positive integer"))
  end
  return n
end
