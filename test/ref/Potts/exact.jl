const Ns = [2,8]
const Qs = [2,3,6]
const J = 1.0
const Ts = collect(0.5:0.5:10.0)

function exact(Q, J, N, Ts)
    j = J/N
    Es = zeros(0)
    spins = zeros(Int, N)
    nst = Q^N
    for ist in 0:(nst-1)
        spins .= 0
        for i in 1:N
            spins[i] = mod(ist,Q)
            ist = div(ist,Q)
        end
        E = 0.0
        for q in 1:Q
            nq = count(s->s==q, spins)
            E -= 0.5j*(nq*(nq-1))
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
    for Q in Qs
        Es = exact(Q,J,N,Ts)
        open(@sprintf("Q_%d__J_%.1f__N_%d.dat", Q, J, N), "w") do io
            for (T,E) in zip(Ts, Es)
                @printf(io, "%.15f %.15f\n", T, E)
            end
        end
    end
end
