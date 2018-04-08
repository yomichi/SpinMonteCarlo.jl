const obsnames = ["Energy", "Energy^2", "Specific Heat",
                  # "|Magnetization|",
                  "Magnetization^2",
                  # "Magnetization^4", "Binder Ratio",
                  "Susceptibility",
                  # "Connected Susceptibility",
                 ]

function exact(Ts::AbstractArray, S::Real, Jz::Real, Jxy::Real, Gamma::Real, L::Integer, workdir::AbstractString="work")
    if isdir(workdir)
        rm(workdir, recursive=true)
    end
    mkdir(workdir)
    res = cd(workdir) do 
        runHPhi(S,Jz,Jxy,Gamma,L)
        Es, Szs, Sz2s = load_output()
        finiteT(Ts, Es, Szs, Sz2s, L)
    end
    open(@sprintf("S_%.1f__Jz_%.1f__Jxy_%.1f__G_%.1f__L_%d.dat", S, Jz, Jxy, Gamma, L),"w") do io
        for i in 1:length(Ts)
            @printf(io, "%.15e", Ts[i])
            for n in obsnames
                @printf(io, " %.15e", res[n][i])
            end
            println(io)
        end
    end
end

function runHPhi(S::Real, Jz::Real, Jxy::Real, Gamma::Real, L::Integer)
    open("stdface.def", "w") do io
        println(io, "Lattice = \"chain lattice\"")
        println(io, "Model = \"SpinGC\"")
        println(io, "Method =\"fulldiag\"")
        println(io, "2S = $(round(Int,2S))")
        println(io, "Jz = $Jz")
        println(io, "Jx = $Jxy")
        println(io, "Jy = $Jxy")
        println(io, "Gamma = $Gamma")
        println(io, "L = $L")
    end
    run(`HPhi -s stdface.def`)
end

function load_output()
    Es = zeros(0)
    Szs = zeros(0)
    Sz2s = zeros(0)
    N = 0
    for line in eachline(joinpath("output", "Eigenvalue.dat"))
        words = split(line)
        E = parse(words[2])
        push!(Es, E)
    end
    nst = length(Es)
    for i in 0:(nst-1)
        Sz  = 0.0
        for line in eachline(joinpath("output", "zvo_cisajs_eigen$i.dat"))
            words = split(line)
            site = parse(Int, words[1])
            spin = parse(Int, words[2])
            v = parse(words[5])
            Sz += v*(spin-0.5)
        end
        push!(Szs, Sz)

        Sz2 = 0.0
        for line in eachline(joinpath("output", "zvo_cisajscktalt_eigen$i.dat"))
            words = split(line)
            site_i = parse(Int, words[1])
            spin_i = parse(Int, words[2])
            if parse(Int, words[4]) != spin_i
                continue
            end
            site_j = parse(Int, words[5])
            spin_j = parse(Int, words[6])
            v = parse(words[9])
            Sz2 += ifelse(spin_i==spin_j, v, -v)
        end
        push!(Sz2s, 0.25Sz2)
    end
    
    return Es, Szs, Sz2s
end

function finiteT(Ts::AbstractArray, Es::AbstractArray, Szs::AbstractArray, Sz2s::AbstractArray, N)
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
        sz = Szs[i]/N
        sz2 = Sz2s[i]/(N*N)
        Zs .+= z
        res["Energy"] .+= E.*z
        res["Energy^2"] .+= (E^2).*z
        res["Magnetization"] .+= sz.*z
        # res["|Magnetization|"] .+= abs(sz).*z
        res["Magnetization^2"] .+= sz2.*z
    end
    res["Magnetization"] ./= Zs
    for name in obsnames
        res[name] ./= Zs
    end
    res["Specific Heat"] = N.*(res["Energy^2"] .- res["Energy"].^2).*betas.^2
    res["Susceptibility"] = N.*(res["Magnetization^2"] .- res["Magnetization"].^2).*betas
    # res["Connected Susceptibility"] = N.*(res["Magnetization^2"] .- res["|Magnetization|"].^2).*betas
    # res["Binder Ratio"] = res["Magnetization^4"] ./ (res["Magnetization^2"].^2)
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
