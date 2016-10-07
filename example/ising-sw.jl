include("../src/SpinMonteCarlo.jl")

using SpinMonteCarlo

const Ls = [8, 16, 24]
const Ts = [2.0, 2.1, 2.2, 2.3, 2.4, 2.5]
const MCS = 8192
const Therm = MCS >> 3


for L in Ls
    lat = square_lattice(L)
    model = Ising(lat)
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
            m2, m4 = magnetizations(swinfo, nsites)
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
