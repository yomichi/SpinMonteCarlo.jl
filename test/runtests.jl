using Compat
using Compat.Test

using SpinMonteCarlo

const SEED = 137
const SEED2 = 19937
const MCS = 10000
const Therm = MCS
const alpha = 0.001

@testset begin
    for filename in ("lattice.jl",
                     "classical.jl",
                     "quantum.jl",
                     "checkpoint.jl",
                    )
        t = @elapsed include(filename)
        println("$(filename): $t sec")
    end
end
