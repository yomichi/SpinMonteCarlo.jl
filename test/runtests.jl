using Base.Test

include("../src/SpinMonteCarlo.jl")
using SpinMonteCarlo

const SEED = 19937
const MCS = 20000
const Therm = MCS
const alpha = 0.01

@testset begin
    for filename in ("lattice.jl",
                     "anisotropic.jl",
                     "classical.jl",
                     "quantum.jl",
                    )
        t = @elapsed include(filename)
        println("$(filename): $t sec")
    end
end
