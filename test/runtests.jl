using Base.Test

include("../src/SpinMonteCarlo.jl")
using SpinMonteCarlo

const SEED = 19937
const MCS = 100000
const Therm = 10000
const alpha = 0.01
const obsnames_ising = ["Magnetization", "|Magnetization|", "Magnetization^2", "Magnetization^4", "Binder Ratio",
                        "Susceptibility", "Connected Susceptibility",
                        "Energy", "Energy^2", "Specific Heat",
                       ]
const obsnames_potts = obsnames_ising

const obsnames_xy = ["|Magnetization|", "|Magnetization|^2", "|Magnetization|^4",
                        "Binder Ratio", "Susceptibility", "Connected Susceptibility",
                        "Magnetization x", "|Magnetization x|", "Magnetization x^2", "Magnetization x^4",
                        "Binder Ratio x", "Susceptibility x", "Connected Susceptibility x",
                        "Magnetization y", "|Magnetization y|", "Magnetization y^2", "Magnetization y^4",
                        "Binder Ratio y", "Susceptibility y", "Connected Susceptibility y",
                        "Helicity Modulus x", "Helicity Modulus y",
                        "Energy", "Energy^2", "Specific Heat",
                       ]
const obsnames_clock = obsnames_xy

@testset begin
    include("lattice.jl")
    include("dimer.jl")
    include("square.jl")
    include("anisotropic.jl")
    include("QuantumXXZ.jl")
end
