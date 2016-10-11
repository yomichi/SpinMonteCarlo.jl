include("../src/SpinMonteCarlo.jl")

using SpinMonteCarlo

const Ls = [16, 24, 32]
const Q = 10
const Tc = 1.0/log1p(sqrt(Q))
const Ts = Tc*linspace(0.8, 1.2, 21)
const MCS = 8192
const Therm = MCS >> 3


for L in Ls
    lat = square_lattice(L)
    model = Potts(lat)
    nsites = numsites(lat)
    for T in Ts
        for i in 1:Therm
            SW_update!(model, T)
        end
        obs = BinningObservableSet()
        makeMCObservable!(obs, "Magnetization^2")
        makeMCObservable!(obs, "Magnetization^4")
        for i in 1:MCS
            swinfo = SW_update!(model, T)
            m2, m4 = magnetizations(swinfo, model)
            obs["Magnetization^2"] << m2
            obs["Magnetization^4"] << m4
        end
        jk = jackknife(obs)
        jk["Binder Ratio"] = jk["Magnetization^4"] / (jk["Magnetization^2"] ^ 2)
        @printf("%d %lf %lf %lf %lf %lf %lf %lf\n",
                L, T, 
                mean(jk["Magnetization^2"]), stderror(jk["Magnetization^2"]),
                mean(jk["Magnetization^4"]), stderror(jk["Magnetization^4"]),
                mean(jk["Binder Ratio"]), stderror(jk["Binder Ratio"]))
    end
    println()
    println()
end
