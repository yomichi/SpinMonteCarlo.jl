using SpinMonteCarlo
using Printf

const model = Ising
const lat = "square lattice"
const L = 16
const update = SW_update!

const Tc = 2.0/log1p(sqrt(2))
const Ts = Tc*range(0.85, stop=1.15, length=31)
const MCS = 8192
const Therm = MCS >> 3

for T in Ts
    param = Parameter("Model"=>model, "Lattice"=>lat,
                      "L"=>L, "T"=>T, "J"=>1.0,
                      "Update Method"=>update,
                      "MCS"=>MCS, "Thermalization"=>Therm,
                     )
    result = runMC(param)
    println(@sprintf("%f %.15f %.15f",
                      T, mean(result["Specific Heat"]), stderror(result["Specific Heat"])))
end
