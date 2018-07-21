const obsnames = ["Energy", "Energy^2", "Specific Heat",
                  "|Magnetization x|", "Magnetization x^2", "Magnetization x^4",
                  "Susceptibility x", "Connected Susceptibility x", "Binder Ratio x",
                  "|Magnetization y|", "Magnetization y^2", "Magnetization y^4",
                  "Susceptibility y", "Connected Susceptibility y", "Binder Ratio y",
                  "|Magnetization|", "|Magnetization|^2", "|Magnetization|^4",
                  "Susceptibility", "Connected Susceptibility", "Binder Ratio",
                 ]
const Ns = [2,3,4]
const J = 1.0
const Ts = [0.3, 1.0, 3.0, 10.0]

"""
exact finite-T energy (per site) of XY model on `N` site fully connected network
"""
function integrate(J, N, Ts; Q::Integer=50)
    j = J/N
    Es = zeros(0)
    Ms = Vector{Float64}[]
    spins = zeros(Int, N)
    nst = Q^N
    for ist in 0:(nst-1)
        spins .= 1
        for i in 1:N
            spins[i] = mod(ist,Q) + 1
            ist = div(ist,Q)
        end
        E = 0.0
        M = zeros(2)
        nqs = [count(s->s==q, spins) for q in 1:Q]
        for q in 1:Q
            nq = nqs[q]
            E -= 0.5j*(nq*(nq-1))
            for q2 in (q+1):Q
                nq2 = nqs[q2]
                E -= j*nq*nq2*cospi(2*(q-q2)/Q)
            end
            M[1] += nq*cospi(2*q/Q)
            M[2] += nq*sinpi(2*q/Q)
        end
        M ./= N
        push!(Es, E)
        push!(Ms, M)
    end
    sorted = sortperm(Es, rev=true)
    E0 = Es[sorted[end]]

    mbetas = Ts .\ (-1.0)
    nT = length(Ts)
    res = Dict(n=>zeros(nT) for n in obsnames)
    res["Magnetization x"] = zeros(nT)
    res["Magnetization y"] = zeros(nT)
    Zs = zeros(nT)
    for i in sorted
        E = Es[i] - E0
        zs = exp.(mbetas.*E)
        Zs .+= zs
        E += E0
        res["Energy"] .+= zs.*(E/N)
        res["Energy^2"] .+= zs.*(E/N)^2

        M = Ms[i]
        for (d,c) in ((1,"x"),(2,"y"))
            m = M[d]
            res["Magnetization $(c)"] .+= zs.*m
            res["|Magnetization $(c)|"] .+= zs.*abs(m)
            res["Magnetization $(c)^2"] .+= zs.*m^2
            res["Magnetization $(c)^4"] .+= zs.*m^4
        end
        m2 = sum(abs2,M)
        res["|Magnetization|"] .+= zs.*sqrt(m2)
        res["|Magnetization|^2"] .+= zs.*m2
        res["|Magnetization|^4"] .+= zs.*m2^2
    end
    res["Magnetization x"] ./= Zs
    res["Magnetization y"] ./= Zs
    for n in obsnames
        res[n] ./= Zs
    end
    return res
end

betas = 1.0 ./ Ts

for N in Ns
    res = integrate(J,N,Ts)
    res["Specific Heat"] = N.*(res["Energy^2"] .- res["Energy"].^2).*betas.^2
    for c in ("x","y")
        res["Susceptibility $c"] = N.*(res["Magnetization $c^2"] .- res["Magnetization $c"].^2).*betas
        res["Connected Susceptibility $c"] = N.*(res["Magnetization $c^2"] .- res["|Magnetization $c|"].^2).*betas
        res["Binder Ratio $c"] = res["Magnetization $c^4"] ./ (res["Magnetization $c^2"].^2)
    end
    res["Susceptibility"] = N.*res["|Magnetization|^2"].*betas
    res["Connected Susceptibility"] = N.*(res["|Magnetization|^2"] .- res["|Magnetization|"].^2).*betas
    res["Binder Ratio"] = res["|Magnetization|^4"] ./ (res["|Magnetization|^2"].^2)
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
