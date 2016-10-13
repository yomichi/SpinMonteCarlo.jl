include("../src/SpinMonteCarlo.jl")

using SpinMonteCarlo

const Ls = [8, 16, 24, 32]
const Ts = linspace(1.5, 0.5, 11)
const MCS = 8192
const Therm = MCS >> 3


for L in Ls
    lat = square_lattice(L)
    model = XY(lat)
    nsites = numsites(lat)
    for T in Ts
        for i in 1:Therm
            SW_update!(model, T)
        end
        obs = BinningObservableSet()
        makeMCObservable!(obs, "Magnetization^2")
        makeMCObservable!(obs, "Magnetization^4")
        for d in 1:dim(lat)
            makeMCObservable!(obs, "Helicity Modulus $d")
        end
        for i in 1:MCS
            SW_update!(model, T)
            M, E, U = measure(model, T)
            m2 = sumabs2(M)
            obs["Magnetization^2"] << m2
            obs["Magnetization^4"] << m2*m2
            for d in 1:dim(lat)
                obs["Helicity Modulus $d"] << U[d]
            end
        end
        jk = jackknife(obs)
        jk["Binder Ratio"] = jk["Magnetization^4"] / (jk["Magnetization^2"] ^ 2)
        print(L, " ", T, " ")
        dump_plot(STDOUT, jk["Magnetization^2"], put_following_space=true)
        dump_plot(STDOUT, jk["Magnetization^4"], put_following_space=true)
        dump_plot(STDOUT, jk["Binder Ratio"], put_following_space=true)
        for d in 1:dim(lat)
            dump_plot(STDOUT, jk["Helicity Modulus $d"], put_following_space=true)
        end
        println()
    end
    println()
    println()
end
