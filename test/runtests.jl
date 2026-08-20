using Test

using SpinMonteCarlo

const SEED = 137
const SEED2 = 19937
const MCS = 10000
const Therm = MCS
const alpha = 0.001

@testset begin
    filenames = ["aqua.jl",
                 "observable.jl",
                 "lattice.jl",
                 "rng.jl",
                 "classical.jl",
                 "quantum.jl",
                 "checkpoint.jl",
                 "snapshot.jl",
                 "union_find.jl"]
    for filename in filenames
        t = @elapsed include(filename)
        println("$(filename): $t sec")
    end
end
