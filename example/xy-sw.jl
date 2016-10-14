include("../src/SpinMonteCarlo.jl")

using SpinMonteCarlo

const Ls = [8, 16, 24, 32]
const Ts = linspace(1.5, 0.5, 11)
const MCS = 8192
const Therm = MCS >> 3

const io = open("res-xy-sw.dat", "w")

for (i, name) in enumerate(["L", "T",
                            "Seconds per step", "Steps per second",
                            "M^2", "error of M^2",
                            "M^4", "error of M^4", 
                            "Binder ratio", "error of Binder ratio", 
                            ])
    println(io, "# \$$i : $name")
end

println(io, "# \$11- : Helicity modulus")


for L in Ls
    lat = square_lattice(L)
    model = XY(lat)
    nsites = numsites(lat)
    for T in Ts
        for i in 1:Therm
            SW_update!(model, T)
        end
        obs = BinningObservableSet()
        makeMCObservable!(obs, "Time")
        makeMCObservable!(obs, "Speed")
        makeMCObservable!(obs, "Magnetization^2")
        makeMCObservable!(obs, "Magnetization^4")
        for d in 1:dim(lat)
            makeMCObservable!(obs, "Helicity Modulus $d")
        end
        for i in 1:MCS
            tic()
            SW_update!(model, T)
            M, E, U = measure(model, T)
            t = toq()

            m2 = sumabs2(M)
            obs["Time"] << t
            obs["Speed"] << 1.0/t
            obs["Magnetization^2"] << m2
            obs["Magnetization^4"] << m2*m2
            for d in 1:dim(lat)
                obs["Helicity Modulus $d"] << U[d]
            end
        end
        jk = jackknife(obs)
        jk["Binder Ratio"] = jk["Magnetization^4"] / (jk["Magnetization^2"] ^ 2)
        print(io, L, " ", T, " ")
        print(io, mean(jk["Time"]), " ", mean(jk["Speed"]), " ")
        dump_plot(io, jk["Magnetization^2"], put_following_space=true)
        dump_plot(io, jk["Magnetization^4"], put_following_space=true)
        dump_plot(io, jk["Binder Ratio"], put_following_space=true)
        for d in 1:dim(lat)
            dump_plot(io, jk["Helicity Modulus $d"], put_following_space=true)
        end
        println(io)
    end
    println(io)
    println(io)
    info("L = $L finished.")
end
close(io)
