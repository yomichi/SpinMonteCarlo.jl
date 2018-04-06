const Ns = [2,8]
const J = 1.0
const Ts = [0.5, 1.0, 2.0, 5.0, 10.0]

"""
Exact finite-T energy (per site) for Ising model on `N` site fully connected network
"""
function exact(J, N, Ts)
    j = J/N
    Es = zeros(N+1)
    ns = zeros(Int,N+1)
    for Nup in 0:N
        Ndown = N - Nup
        n = binomial(N, Nup)
        push!(ns, n)
        E = j*(Nup*Ndown - 0.5(Nup*(Nup-1)) - 0.5(Ndown*(Ndown-1)))
        push!(Es, E)
    end
    sorted = sortperm(Es, rev=true)
    E0 = Es[sorted[end]]

    mbetas = Ts .\ (-1.0)
    nT = length(Ts)
    Zs = zeros(nT)
    EZs = zeros(nT)
    for i in sorted
        E = Es[i] - E0
        n = ns[i]
        zs = exp.(mbetas.*E) .* n
        Zs .+= zs
        EZs .+= zs.*(E/N)
    end
    return EZs ./ Zs .+ (E0/N)
end

for N in Ns
    Es = exact(J,N,Ts)
    open(@sprintf("J_%.1f__N_%d.dat", J, N), "w") do io
        for (T,E) in zip(Ts, Es)
            @printf(io, "%.15f %.15f\n", T, E)
        end
    end
end
