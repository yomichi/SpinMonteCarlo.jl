include("../src/SpinMonteCarlo.jl")

using SpinMonteCarlo

const J = -1.0
const Ls = [8, 16, 32]
const gs = linspace(0.0, 1.0, 11)
const MCS = 8192
const Therm = MCS >> 3

params = Dict{String, Any}[]
for L in Ls
    for g in gs
        T = 1.0/L
        push!(params, Dict("Model"=>TransverseFieldIsing,
                           "Lattice"=>chain_lattice,
                           "L"=>L, "T"=>T, "J"=>J, "Gamma"=>g,
                           "MCS"=>MCS, "Therm"=>Therm, "Verbose"=>true,
                          ))
    end
end

obs = map(runMC, params)

const pnames = ["L", "T", "Gamma"]
const onames = ["Magnetization",
                "|Magnetization|",
                "Magnetization^2",
                "Magnetization^4",
                "Binder Ratio",
                "Susceptibility",
                "Connected Susceptibility",
                "MCS per Second",
                "Time per MCS",
               ]

const io = open("res-TFI.dat", "w")
i=1
for pname in pnames
    println(io, "# \$$i : $pname")
    i+=1
end
for oname in onames
    println(io, "# \$$i, $(i+1): $oname")
    i+=2
end

for (p,o) in zip(params, obs)
    for pname in pnames
        print(io, p[pname], " ")
    end
    for oname in onames
        @printf(io, "%.15f %.15f ", mean(o[oname]), stderror(o[oname]))
    end
    println(io)
end

