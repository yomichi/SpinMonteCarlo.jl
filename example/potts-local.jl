include("../src/SpinMonteCarlo.jl")

using SpinMonteCarlo

const Ls = [8, 16, 24]
const Qs = [3, 10]
const MCS = 8192
const Therm = MCS >> 3

const io = open("res-potts-local.dat", "w")

for (i, name) in enumerate(["Q", "L", "T",
                            "Seconds per step", "Steps per second",
                            "M^2", "error of M^2",
                            "M^4", "error of M^4", 
                            "Binder ratio", "error of Binder ratio", 
                            ])
    println(io, "# \$$i : $name")
end

for Q in Qs
    Tc = 1.0/log1p(sqrt(Q))
    Ts = Tc*linspace(0.8, 1.2, 21)
    for L in Ls
        lat = square_lattice(L)
        model = Potts(lat, Q)
        nsites = numsites(lat)
        for T in Ts
            for i in 1:Therm
                local_update!(model, T)
            end
            obs = BinningObservableSet()
            makeMCObservable!(obs, "Time")
            makeMCObservable!(obs, "Speed")
            makeMCObservable!(obs, "Magnetization^2")
            makeMCObservable!(obs, "Magnetization^4")
            for i in 1:MCS
                tic()
                local_update!(model, T)
                M, E = measure(model, T)
                m2 = M*M
                m4 = m2*m2
                t = toq()
                obs["Time"] << t
                obs["Speed"] << 1.0/t
                obs["Magnetization^2"] << m2
                obs["Magnetization^4"] << m4
            end
            jk = jackknife(obs)
            jk["Binder Ratio"] = jk["Magnetization^4"] / (jk["Magnetization^2"] ^ 2)

            print(io, Q, " ", L, " ", T, " ")
            print(io, mean(jk["Time"]), " ", mean(jk["Speed"]), " ")
            dump_plot(io, jk["Magnetization^2"], put_following_space=true)
            dump_plot(io, jk["Magnetization^4"], put_following_space=true)
            dump_plot(io, jk["Binder Ratio"], put_following_space=true)
            println(io)
        end
        println(io)
        println(io)
        info("L = $L finished.")
    end
    info("Q = $Q finished.")
end
close(io)