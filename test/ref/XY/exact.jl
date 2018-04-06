const Ns = [2,3,4]
const J = 1.0
const Ts = [0.2, 0.5, 1.0, 2.0, 5.0, 10.0]
# const Ts = linspace(0.01, 1.00, 100)

"""
exact finite-T energy (per site) of XY model on `N` site fully connected network
"""
function integrate(J, N, Ts; Q::Integer=50)
    j = J/N
    Es = zeros(0)
    spins = zeros(Int, N)
    nst = Q^N
    for ist in 0:(nst-1)
        if ist%10000==0
            @show ist/nst
        end
        spins .= 1
        for i in 1:N
            spins[i] = mod(ist,Q) + 1
            ist = div(ist,Q)
        end
        E = 0.0
        nqs = [count(s->s==q, spins) for q in 1:Q]
        for q in 1:Q
            nq = nqs[q]
            E -= 0.5j*(nq*(nq-1))
            for q2 in (q+1):Q
                nq2 = nqs[q2]
                E -= j*nq*nq2*cospi(2*(q-q2)/Q)
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
    Es = integrate(J,N,Ts)
    open(@sprintf("J_%.1f__N_%d.dat", J, N), "w") do io
        for (T,E) in zip(Ts, Es)
            @printf(io, "%.15f %.15f\n", T, E)
        end
    end
end
