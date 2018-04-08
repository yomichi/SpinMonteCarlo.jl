const obsnames = ["Energy", "Energy^2", "Specific Heat",
                  "Magnetization^2", "Magnetization^4", 
                  "Binder Ratio", "Susceptibility",
                 ]

isvalidspin(S::Real) = round(2S)==2S

function genSz(S::Real=0.5)
    if !isvalidspin(S)
        error("invalid spin size: S = $S")
    end
    S2 = Int(2S)
    return diagm([0.5s for s in S2:-2:(-S2)])
end

function genSplus(S::Real=0.5)
    if !isvalidspin(S)
        error("invalid spin size: S = $S")
    end
    sz = genSz(S)
    N = Int(2S+1)
    Splus = zeros(N,N)
    for i in 1:(N-1)
        j = i+1
        m = sz[j,j]
        Splus[i,j] = sqrt((S-m)*(S+m+1))
    end
    return Splus
end

function genSminus(S::Real=1)
    if !isvalidspin(S)
        error("invalid spin size: S = $S")
    end
    sz = genSz(S)

    N = Int(2S+1)
    Sminus = zeros(N,N)
    for i in 2:N
        j = i-1
        m = sz[j,j]
        Sminus[i,j] = sqrt((S+m)*(S-m+1))
    end
    return Sminus
end

genSx(S::Real=0.5) = 0.5.*(genSplus(S) .+ genSminus(S))

function genSztotal(S::Real, L::Integer)
    if !isvalidspin(S)
        error("invalid spin size: S = $S")
    end
    S2 = Int(2S)
    Nlocal = S2+1
    sz = genSz(S)

    res = zeros(Nlocal^L, Nlocal^L)
    for i in 1:L
        res .+= kron(kron(eye(Nlocal^(i-1)), sz), eye(Nlocal^(L-i)))
    end
    return res
end

function genH(S::Real, Jz::Real, Jxy::Real, Gamma::Real, L::Integer)
    if !isvalidspin(S)
        error("invalid spin size: S = $S")
    end
    S2 = Int(2S)
    Nlocal = S2+1

    sz = genSz(S)
    sx = genSx(S)
    sp = genSplus(S)
    sm = genSminus(S)
    ss = Jz.*kron(sz,sz) .+ (0.5Jxy).*(kron(sp,sm) .+ kron(sm,sp)) 

    res = Jz.*kron(kron(sz, eye(Nlocal^(L-2))), sz)
    res .+= (0.5Jxy).*kron(kron(sm, eye(Nlocal^(L-2))), sp)
    res .+= (0.5Jxy).*kron(kron(sp, eye(Nlocal^(L-2))), sm)

    for i in 1:(L-1)
        res .+= kron(kron(eye(Nlocal^(i-1)), ss), eye(Nlocal^(L-i-1)))
    end

    for i in 1:L
        res .-= Gamma .* kron(kron(eye(Nlocal^(i-1)), sx), eye(Nlocal^(L-i)))
    end

    return res
end

function exact(Ts::AbstractArray, S::Real, Jz::Real, Jxy::Real, Gamma::Real, L::Integer)
    Es, Ms, M2s, M4s = diagonalize(S,Jz,Jxy,Gamma,L)
    res = finiteT(Ts, Es, Ms, M2s, M4s, L)
    filename = @sprintf("S_%.1f__Jz_%.1f__Jxy_%.1f__G_%.1f__L_%d.dat", S, Jz, Jxy, Gamma, L)
    open(filename,"w") do io
        for i in 1:length(Ts)
            @printf(io, "%.15e", Ts[i])
            for n in obsnames
                @printf(io, " %.15e", res[n][i])
            end
            println(io)
        end
    end
    info("$filename done")
end

function diagonalize(S::Real, Jz::Real, Jxy::Real, Gamma::Real, L::Integer)
    H = genH(S, Jz, Jxy, Gamma, L)
    sz = genSztotal(S, L)
    Es, vecs = eig(H)
    V = sz*vecs
    Ms = diag(vecs'*V)
    V .= sz*V
    M2s = diag(vecs'*V)
    V .= sz*sz*V
    M4s = diag(vecs'*V)
    return Es, Ms./L, M2s./(L^2), M4s./(L^4)
end

function finiteT(Ts::AbstractArray, Es::AbstractArray, Ms::AbstractArray, M2s::AbstractArray, M4s::AbstractArray, N)
    nT = length(Ts)
    betas = 1.0./Ts
    Zs = zeros(nT)
    res = Dict(name=>zeros(nT) for name in obsnames)
    res["Magnetization"] = zeros(nT)
    sorted = sortperm(Es, rev=true)
    E0 = Es[sorted[end]]
    for i in sorted
        E = Es[i] - E0
        z = exp.(betas.*(-E))
        E += E0
        E /= N
        Zs .+= z
        res["Energy"] .+= E.*z
        res["Energy^2"] .+= (E^2).*z
        res["Magnetization"] .+= Ms[i].*z
        res["Magnetization^2"] .+= M2s[i]*z
        res["Magnetization^4"] .+= M4s[i]*z
    end
    res["Magnetization"] ./= Zs
    for name in obsnames
        res[name] ./= Zs
    end
    res["Specific Heat"] = N.*(res["Energy^2"] .- res["Energy"].^2).*betas.^2
    res["Susceptibility"] = N.*(res["Magnetization^2"] .- res["Magnetization"].^2).*betas
    res["Binder Ratio"] = res["Magnetization^4"] ./ (res["Magnetization^2"].^2)
    return res
end

const Ts = collect(1.0:1.0:10.01)
exact(Ts, 0.5, -1.0, 0.0, 0.0, 8)
exact(Ts, 0.5, -1.0, 0.5, 0.0, 8)
exact(Ts, 0.5, -1.0, 1.0, 0.0, 8)
exact(Ts, 0.5, -1.0, 2.0, 0.0, 8)
exact(Ts, 0.5,  1.0, 0.0, 0.0, 8)
exact(Ts, 0.5,  1.0, 0.5, 0.0, 8)
exact(Ts, 0.5,  1.0, 1.0, 0.0, 8)
exact(Ts, 0.5,  1.0, 2.0, 0.0, 8)
exact(Ts, 0.5,  0.0, 1.0, 0.0, 8)
exact(Ts, 0.5, -1.0, 0.0, 1.0, 8)
exact(Ts, 0.5,  0.0, 1.0, 0.0, 3)
exact(Ts, 1.0,   1.0, 1.0, 0.0, 6)
exact(Ts, 1.0,   0.0, 1.0, 0.0, 6)
exact(Ts, 1.0,  -1.0, 0.0, 1.0, 6)
