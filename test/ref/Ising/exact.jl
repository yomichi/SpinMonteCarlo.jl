const Ns = [2,8]
const J = 1.0
const Ts = collect(0.5:0.5:10.0)

function exact(J, N, Ts)
    j = J/N
    nT = length(Ts)
    mbetas = Ts .\ (-1.0)
    Zs = zeros(nT)
    EZs = zeros(nT)
    for Nup in 0:N
        Ndown = N - Nup
        n = binomial(N, Nup)
        E = j*(Nup*Ndown - 0.5(Nup*(Nup-1)) - 0.5(Ndown*(Ndown-1)))
        zs = exp.(mbetas.*E) .* n
        Zs .+= zs
        EZs .+= zs.*(E/N)
    end
    return EZs ./ Zs
end

for N in Ns
    Es = exact(J,N,Ts)
    open(@sprintf("J_%.1f__N_%d.dat", J, N), "w") do io
        for (T,E) in zip(Ts, Es)
            @printf(io, "%.15f %.15f\n", T, E)
        end
    end
end
