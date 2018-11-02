using SpinMonteCarlo
using Printf

const ups = [0,1,2]
const Qs = [2,3,4]
const Ls = [8,16,24]
const MCS = 8192
const Therm = MCS >> 3

params = Dict{String, Any}[]
for update in ups
    for Q in Qs
        const Tc = 1.0/log1p(sqrt(Q))
        const Ts = Tc*range(0.85, stop=1.15, length=11)
        for L in Ls
            for T in Ts
                push!(params,
                  Parameter("Model"=>Potts, "Lattice"=>"square lattice",
                            "Q"=>Q, "L"=>L, "T"=>T, "J"=>1.0,
                            "MCS"=>MCS, "Thermalization"=>Therm,
                            "Update Method"=> (update==0 ? local_update! :
                                               update==1 ? SW_update! : Wolff_update!),
                            "update"=>update,
                            "Verbose"=>true,
                           ))
            end
        end
    end
end

obs = map(runMC, params)

const pnames = ["update", "Q", "L", "T"]
const onames = ["Magnetization",
                "|Magnetization|",
                "Magnetization^2",
                "Magnetization^4",
                "Binder Ratio",
                "Susceptibility",
                "Connected Susceptibility",
                "Energy",
                "Energy^2",
                "Specific Heat",
                "MCS per Second",
                "Time per MCS",
               ]

const io = open("res-potts.dat", "w")
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

