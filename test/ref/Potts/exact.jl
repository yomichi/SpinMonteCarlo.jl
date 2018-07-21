const obsnames = ["Energy", "Energy^2", "Specific Heat",
                  "|Magnetization|", "Magnetization^2", "Magnetization^4",
                  "Susceptibility", "Connected Susceptibility", "Binder Ratio",
                 ]
const Ns = [2,8]
const Qs = [2,4]
const J = 1.0
const Ts = [0.3, 1.0, 3.0, 10.0]

"""
exact calculation for `Q` state Potts model on `N` site fully connected network
"""
function exact(Q, J, N, Ts)
    j = J/N
    Es = zeros(0)
    Ms = zeros(0)
    spins = zeros(Int, N)
    nst = Q^N
    for ist in 0:(nst-1)
        spins .= 1
        for i in 1:N
            spins[i] = mod(ist,Q)+1
            ist = div(ist,Q)
        end
        E = 0.0
        M = 0.0
        Qinv = 1.0/Q
        nq = count(s->s==1, spins)
        E -= 0.5j*(nq*(nq-1))
        M += nq*(1.0-Qinv)
        for q in 2:Q
            nq = count(s->s==q, spins)
            E -= 0.5j*(nq*(nq-1))
            M -= nq*Qinv
        end
        M /= N
        push!(Es, E)
        push!(Ms, M)
    end
    sorted = sortperm(Es,rev=true)
    E0 = Es[sorted[end]]

    mbetas = Ts .\ (-1.0)
    nT = length(Ts)
    res = Dict(n=>zeros(nT) for n in obsnames)
    res["Magnetization"] = zeros(nT)
    Zs = zeros(nT)
    for i in sorted
        E = Es[i] - E0
        zs = exp.(mbetas.*E)
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
    for Q in Qs
        res = exact(Q,J,N,Ts)
        res["Specific Heat"] = N.*(res["Energy^2"] .- res["Energy"].^2).*betas.^2
        res["Susceptibility"] = N.*(res["Magnetization^2"] .- res["Magnetization"].^2).*betas
        res["Connected Susceptibility"] = N.*(res["Magnetization^2"] .- res["|Magnetization|"].^2).*betas
        res["Binder Ratio"] = res["Magnetization^4"] ./ (res["Magnetization^2"].^2)
        open(@sprintf("Q_%d__J_%.1f__N_%d.dat", Q, J, N), "w") do io
            for i in 1:length(Ts)
                @printf(io, "%.15f", Ts[i])
                for n in obsnames
                    @printf(io, " %.15f", res[n][i])
                end
                println(io)
            end
        end
    end
end
