using Base.Test

include("../src/SpinMonteCarlo.jl")
using SpinMonteCarlo

const SEED = 19937

include("dimer.jl")
include("updatemethod.jl")
