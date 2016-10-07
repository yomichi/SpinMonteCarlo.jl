using LsqFit

import Base: mean, var, sum

export BinningVectorObservable, push!, tau, reset!
export extrapolate_tau, extrapolate_stderror
export BinningVectorObservableSet

type BinningVectorObservable <: VectorObservable
  ## index of these vectors denotes the level of bins (each bin stores the mean of 2^(i-1) values)
  raw_ts :: Vector{Vector{Float64}}        ## Time series of raw data
  bins :: Vector{Vector{Float64}}          ## Time series of bins
  sum :: Vector{Vector{Float64}}           ## Summation of stored values
  sum2 :: Vector{Vector{Float64}}          ## Summation of square of bin's mean
  entries :: Vector{Int}           ## Number of bins
  binsize :: Int
  lastbin :: Int
  minbinnum :: Int
  maxlevel :: Int
end

BinningVectorObservable(minbinnum::Int = 128) = BinningVectorObservable( 
  Vector{Float64}[], Vector{Float64}[], Vector{Float64}[], Vector{Float64}[], zeros(Int,1), 1, 0, minbinnum, 1
)

function reset!(b::BinningVectorObservable)
  b.raw_ts = Vector{Float64}[]
  b.bins = Vector{Float64}[]
  b.sum = Vector{Float64}[]
  b.sum2 = Vector{Float64}[]
  b.entries = zeros(Int,1)
  b.binsize = 1
  b.lastbin = 0
  b.maxlevel = 1

  return b
end

maxlevel(b::BinningVectorObservable) = b.maxlevel
count(b::BinningVectorObservable, level::Int = 1) = b.entries[level]
sum(b::BinningVectorObservable, level::Int = 1) = b.sum[level]
sum2(b::BinningVectorObservable, level::Int = 1) = b.sum2[level]

function push!(b::BinningVectorObservable, x::Vector)
  ## time series
  push!(b.raw_ts, deepcopy(x))
  if b.lastbin == 0
    push!(b.bins, deepcopy(x))
    b.lastbin += 1
  else
    b.bins[end] += x
    b.lastbin += 1
  end
  if b.lastbin == b.binsize
    b.bins[end] /= b.binsize
    b.lastbin = 0
    if length(b.bins) == 2(b.minbinnum)
      new_bins = Vector{Vector{Float64}}(b.minbinnum)
      for i in 1:b.minbinnum
        new_bins[i] = 0.5(b.bins[2i-1]+b.bins[2i])
      end
      b.bins = new_bins
      b.binsize <<= 1
      b.maxlevel += 1
    end
  end

  ## binning
  if b.entries[1] == 0
    push!(b.sum, deepcopy(x))
    push!(b.sum2, squared(x))
  else
    b.sum[1] += x
    b.sum2[1] += squared(x)
  end
  i = b.entries[1]
  b.entries[1] += 1
  level = 2
  bsize = 2
  while i & 1 > 0
    if level > length(b.sum)
      push!(b.sum, zeros(length(x)))
      push!(b.sum2,zeros(length(x)))
      push!(b.entries, 0)
    end
    lastbin = b.sum[1]/bsize - b.sum[level]
    if b.entries[level] == 0
      b.sum[level] = lastbin
      b.sum2[level] = squared(lastbin)
      b.entries[level] += 1
    else
      b.sum[level] += lastbin
      b.sum2[level] += squared(lastbin)
      b.entries[level] += 1
    end

    level += 1
    bsize <<= 1
    i >>= 1
  end

  return b
end

function mean(b::BinningVectorObservable, level::Int = 1)
  return sum(b, level)/count(b, level)
end

function var(b::BinningVectorObservable, level::Int = 1)
  n = count(b, level)
  s = sum(b, level)
  s2 = sum2(b, level)
  if n > 1
    v2 = s2 - s.*s/n
    return map(maxzero, v2/(n-1))
  elseif n == 1
    return fill(Inf, length(b.sum[1]))
  else
    return NaN
  end
end

function stderror(b::BinningVectorObservable, level::Int = maxlevel(b))
  return sqrt(var(b,level)/count(b,level))
end

function confidence_interval(b::BinningVectorObservable, confidence_rate::Real, level::Int = maxlevel(b))
  q = 0.5+0.5*confidence_rate
  correction = quantile( TDist(count(b,level)), q)
  serr = stderror(b, level)
  return correction * serr
end
function confidence_interval(b::BinningVectorObservable, confidence_rate_symbol::Symbol = :sigma1, level::Int = maxlevel(b))
  n = parsesigma(confidence_rate_symbol)
  return confidence_interval(b, erf(0.5n*sqrt(2.0)), level)
end

function tau(b::BinningVectorObservable, level::Int = maxlevel(b))
  binsize = 1<<(level-1)
  return 0.5*( (binsize*var(b,level))./var(b) - 1.0)
end

#=
linearmodel(x::Float64, p::Vector{Float64}) = p[1] + x*p[2]
linearmodel(xs::Vector{Float64}, p::Vector{Float64}) = map(x->linearmodel(x,p),xs)
=#

function extrapolate_detail(op :: Function, b::BinningVectorObservable, point::Int)
  ml = maxlevel(b)
  ll = max(ml-point+1, 1)
  levels = ll:ml
  ns = map( level->1<<(level-1), levels)
  ninvs = ns .\ 1.0
  ys = map( level->op(b, level), levels)
  fit = curve_fit(linearmodel, ninvs, ys, [ys[end], 0.0])
  return fit.param[1], estimate_errors(fit)[1]
end
extrapolate_tau(b::BinningVectorObservable, point::Int = 5) = extrapolate_detail(tau, b, point)
extrapolate_stderror(b::BinningVectorObservable, point::Int = 5) = extrapolate_detail(stderror, b, point)

function show(io::IO, obs::BinningVectorObservable)
  if count(obs) > 0
    print(io, mean(obs), " +/- ", stderror(obs), "; tau = ", tau(obs))
  else
    print(io, "No entries")
  end
end


typealias BinningVectorObservableSet MCObservableSet{BinningVectorObservable}

