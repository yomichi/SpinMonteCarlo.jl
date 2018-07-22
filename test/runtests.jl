using Compat
using Compat.Test

using SpinMonteCarlo

const SEED = 137
const SEED2 = 19937
const MCS = 10000
const Therm = MCS
const alpha = 0.001

@testset begin
    filenames = ["lattice.jl",
                 "classical.jl",
                 "quantum.jl",
                ]
    if VERSION <= v"0.6.4"
        push!(filenames, "checkpoint.jl")
    end
    for filename in filenames
        t = @elapsed include(filename)
        println("$(filename): $t sec")
    end
end
