function loaddata(filename)
    Ts = zeros(0)
    Es = zeros(0)
    for line in eachline(joinpath("ref", filename))
        words = split(line)
        push!(Ts, parse(words[1]))
        push!(Es, parse(words[2]))
    end
    return Ts, Es
end

function parse_filename(filename)
    m = match(r"^S_([\d.-]*)__Jz_([\d.-]*)__Jxy_([\d.-]*)__G_([\d.-]*)__L_([\d.-]*).dat$", filename)
    if m == nothing
        return nothing
    end
    p = Dict()
    p[:S] = parse(m.captures[1])
    p[:Jz] = parse(m.captures[2])
    p[:Jxy] = parse(m.captures[3])
    p[:Gamma] = parse(m.captures[4])
    p[:L] = parse(Int, m.captures[5])
    return p
end

function QMC(T; S=0.5, Jz=1.0, Jxy=1.0, Gamma=0.0, L=8)
    p = Dict("Model"=>QuantumXXZ, "Lattice"=>chain_lattice,
             "S"=>S, "L"=>L,
             "Jz"=>Jz, "Jxy"=>Jxy,
             "Gamma"=>Gamma,
             "T"=>T,
             "MCS"=>MCS,
             "Thermalization"=>Therm,
            )
    return runMC(p)
end

@testset "QuantumXXZ chain" begin
    srand(SEED)
    for filename in readdir("ref")
        p = parse_filename(filename)
        if p == nothing
            continue
        end
        @testset "S=$(p[:S]), Jz=$(p[:Jz]), Jxy=$(p[:Jxy]), Gamma=$(p[:Gamma]), L=$(p[:L])" begin
            Ts, exacts = loaddata(filename)
            N = length(Ts)
            for (T,exact) in zip(Ts,exacts)
                res = QMC(T; p...)
                ene = res["Energy"]
                sgn = mean(res["Sign"])
                if !(p_value(ene, exact) > alpha/N)
                    if  sgn < 1.0 && !(isfinite(mean(ene)))
                        continue
                    else
                        @show T
                        @show exact
                        @show ene
                    end
                end
                @test p_value(ene, exact) > alpha/N
            end
        end
    end
end
