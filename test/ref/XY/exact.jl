const Ns = [2,8]
const Qs = [2,6]
const J = 1.0
const Ts = [0.2, 0.5, 1.0, 2.0, 0.5, 10.0]

"""
exact finite-T energy (per site) of XY model on `N` site fully connected network
"""
function exact(J, N, Ts; nsamples::Integer=1_000_000)
    j = J/N
    Es = zeros(0)
    spins = zeros(N)

    for i in 1:nsamples
        spins .= rand(N)
        E = 0.0
        for i in 1:N
            for j in (i+1):N
                E -= j*cospi(2(spins[i]-spins[j]))
            end
        end
        push!(Es, E)
    end
    sort!(Es,rev=true)
    E0 = Es[end]
    Es .-= E0

    mbetas = Ts .\ (-1.0)
    nT = length(Ts)
    Zs = zeros(nT)
    EZs = zeros(nT)
    for E in Es
        z = exp.(mbetas.*E)
        Zs .+= z
        EZs .+= z.*(E/N)
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
