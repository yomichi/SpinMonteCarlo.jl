const obsnames = ["Energy", "Energy^2", "Specific Heat",
                  "|Magnetization|", "Magnetization^2", "Magnetization^4",
                  "Susceptibility", "Connected Susceptibility", "Binder Ratio",
                 ]
const Ns = [2,8]
const J = 1.0
const Ts = [0.3, 1.0, 3.0, 10.0]

"""
Exact calculation for Ising model on `N` site fully connected network
"""
function exact(J, N, Ts)
    j = J/N
    Es = zeros(N+1)
    Ms = zeros(N+1)
    ns = zeros(Int,N+1)
    for Nup in 0:N
        Ndown = N - Nup
        n = binomial(N, Nup)
        push!(ns, n)
        E = j*(Nup*Ndown - 0.5(Nup*(Nup-1)) - 0.5(Ndown*(Ndown-1)))
        push!(Es, E)
        M = (Nup-Ndown)/N
        push!(Ms, M)
    end
    sorted = sortperm(Es, rev=true)
    E0 = Es[sorted[end]]

    mbetas = Ts .\ (-1.0)
    nT = length(Ts)
    res = Dict(n=>zeros(nT) for n in obsnames)
    res["Magnetization"] = zeros(nT)
    Zs = zeros(nT)
    for i in sorted
        E = Es[i] - E0
        zs = exp.(mbetas.*E) * ns[i]
        Zs .+= zs
        E += E0
        res["Energy"] .+= zs.*(E/N)
        res["Energy^2"] .+= zs.*(E/N)^2

        M = Ms[i]
        res["Magnetization"] .+= zs.*M
        res["|Magnetization|"] .+= zs.*abs(M)
        res["Magnetization^2"] .+= zs.*M^2
        res["Magnetization^4"] .+= zs.*M^4
    end
    res["Magnetization"] ./= Zs
    for n in obsnames
        res[n] ./= Zs
    end
    return res
end

betas = 1.0 ./ Ts

for N in Ns
    res = exact(J,N,Ts)
    res["Specific Heat"] = N.*(res["Energy^2"] .- res["Energy"].^2).*betas.^2
    res["Susceptibility"] = N.*(res["Magnetization^2"] .- res["Magnetization"].^2).*betas
    res["Connected Susceptibility"] = N.*(res["Magnetization^2"] .- res["|Magnetization|"].^2).*betas
    res["Binder Ratio"] = res["Magnetization^4"] ./ (res["Magnetization^2"].^2)
    open(@sprintf("J_%.1f__N_%d.dat", J, N), "w") do io
        for i in 1:length(Ts)
            @printf(io, "%.15f", Ts[i])
            for n in obsnames
                @printf(io, " %.15f", res[n][i])
            end
            println(io)
        end
    end
end
