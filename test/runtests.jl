using Base.Test

include("../src/SpinMonteCarlo.jl")
using SpinMonteCarlo

const SEED = 19937

include("lattice.jl")
include("dimer.jl")
include("square.jl")
