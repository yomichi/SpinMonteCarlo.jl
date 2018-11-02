using SpinMonteCarlo
using Printf

const S = 0.5
const J = 1.0
const Gamma = 0.0
const Ls = [8]
const Ts = 0.5:0.5:10.01
const MCS = 8192
const Therm = MCS >> 3

params = Dict{String, Any}[]
for L in Ls
    for T in Ts
        push!(params,
              Parameter("Model"=>QuantumXXZ,
                        "Lattice"=> "chain lattice",
                        "L"=>L, "T"=>T, "J"=>J, "S"=>S, "Gamma"=>Gamma,
                        "Update Method"=>loop_update!,
                        "MCS"=>MCS, "Therm"=>Therm,
                        "Verbose"=>true,
                       ))
    end
end

obs = runMC(params)

const pnames = ["S", "J", "Gamma", "L", "T"]
const onames = ["Magnetization",
                "|Magnetization|",
                "Magnetization^2",
                "Magnetization^4",
                "Binder Ratio",
                "Susceptibility",
                "Connected Susceptibility",
                "Energy",
                "MCS per Second",
                "Time per MCS",
               ]

const io = open("res-QuantumXXZ.dat", "w")
i=1
for pname in pnames
    println(io, "# \$$i : $pname")
    global i+=1
end
for oname in onames
    println(io, "# \$$i, $(i+1): $oname")
    global i+=2
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

