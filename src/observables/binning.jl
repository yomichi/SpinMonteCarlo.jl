using LsqFit

export BinningObservable, tau
export extrapolate_tau, extrapolate_stderror
export BinningObservableSet

type BinningObservable <: ScalarObservable
  ## index of these vectors denotes the level of bins (each bin stores the mean of 2^(i-1) values)
  raw_ts :: Vector{Float64}        ## Time series of raw data
  bins :: Vector{Float64}          ## Time series of bins
  sum :: Vector{Float64}           ## Summation of stored values
  sum2 :: Vector{Float64}          ## Summation of square of bin's mean
  entries :: Vector{Int}           ## Number of bins
  binsize :: Int
  lastbin :: Int
  minbinnum :: Int
  maxlevel :: Int
end

BinningObservable(minbinnum::Int = 128) = BinningObservable( 
  zeros(0), zeros(0), zeros(1), zeros(1), zeros(Int,1), 1, 0, minbinnum, 1
)

function reset!(b::BinningObservable)
  b.raw_ts = zeros(0)
  b.bins = zeros(0)
  b.sum = zeros(1)
  b.sum2 = zeros(1)
  b.entries = zeros(Int,1)
  b.binsize = 1
  b.lastbin = 0
  b.maxlevel = 1

  return b
end

maxlevel(b::BinningObservable) = b.maxlevel
count(b::BinningObservable, level::Int = 1) = b.entries[level]
sum(b::BinningObservable, level::Int = 1) = b.sum[level]
sum2(b::BinningObservable, level::Int = 1) = b.sum2[level]

function push!(b::BinningObservable, x::Real)
  ## time series
  push!(b.raw_ts, x)
  if b.lastbin == 0
    push!(b.bins, x)
    b.lastbin += 1
  else
    b.bins[end] += x
    b.lastbin += 1
  end
  if b.lastbin == b.binsize
    b.bins[end] /= b.binsize
    b.lastbin = 0
    if length(b.bins) == 2(b.minbinnum)
      new_bins = zeros(b.minbinnum)
      for i in 1:b.minbinnum
        new_bins[i] = 0.5(b.bins[2i-1]+b.bins[2i])
      end
      b.bins = new_bins
      b.binsize <<= 1
      b.maxlevel += 1
    end
  end

  ## binning
  b.sum[1] += x
  b.sum2[1] += x*x
  i = b.entries[1]
  b.entries[1] += 1
  level = 2
  bsize = 2
  while i & 1 > 0
    if level > length(b.sum)
      push!(b.sum, 0)
      push!(b.sum2, 0)
      push!(b.entries, 0)
    end
    lastbin = b.sum[1]/bsize - b.sum[level]
    b.sum[level] += lastbin
    b.sum2[level] += lastbin*lastbin
    b.entries[level] += 1

    level += 1
    bsize <<= 1
    i >>= 1
  end

  return b
end

function mean(b::BinningObservable, level::Int = 1)
  return sum(b, level)/count(b, level)
end

function var(b::BinningObservable, level::Int = 1)
  n = count(b, level)
  s = sum(b, level)
  s2 = sum2(b, level)
  if n > 1
    v2 = s2 - s*s/n
    v2 = maxzero(v2)
    return v2/(n-1)
  elseif n == 1
    return inf(Float64)
  else
    return nan(Float64)
  end
end

function stderror(b::BinningObservable, level::Int = maxlevel(b))
  return sqrt(var(b,level)/count(b,level))
end

function confidence_interval(b::BinningObservable, confidence_rate::Real, level::Int = maxlevel(b))
  q = 0.5+0.5*confidence_rate
  correction = quantile( TDist(count(b,level)), q)
  serr = stderror(b, level)
  return correction * serr
end
function confidence_interval(b::BinningObservable, confidence_rate_symbol::Symbol = :sigma1, level::Int = maxlevel(b))
  n = parsesigma(confidence_rate_symbol)
  return confidence_interval(b, erf(0.5n*sqrt(2.0)), level)
end

function tau(b::BinningObservable, level::Int = maxlevel(b))
  binsize = 1<<(level-1)
  return 0.5*( (binsize*var(b,level))/var(b) - 1.0)
end

linearmodel(x::Float64, p::Vector{Float64}) = p[1] + x*p[2]
linearmodel(xs::Vector{Float64}, p::Vector{Float64}) = map(x->linearmodel(x,p),xs)

function extrapolate_detail(op :: Function, b::BinningObservable, point::Int)
  ml = maxlevel(b)
  ll = max(ml-point+1, 1)
  levels = ll:ml
  ns = map( level->1<<(level-1), levels)
  ninvs = ns .\ 1.0
  ys = map( level->op(b, level), levels)
  fit = curve_fit(linearmodel, ninvs, ys, [ys[end], 0.0])
  return fit.param[1], estimate_errors(fit)[1]
end
extrapolate_tau(b::BinningObservable, point::Int = 5) = extrapolate_detail(tau, b, point)
extrapolate_stderror(b::BinningObservable, point::Int = 5) = extrapolate_detail(stderror, b, point)

function show(io::IO, obs::BinningObservable)
  if count(obs) > 0
    print(io, mean(obs), " +/- ", stderror(obs), "; tau = ", tau(obs))
  else
    print(io, "No entries")
  end
end


typealias BinningObservableSet MCObservableSet{BinningObservable}

