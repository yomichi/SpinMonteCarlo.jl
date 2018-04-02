function loaddata(;S2=1, Jz=1.0, Jxy=1.0, G=0.0)
    Ts = zeros(0)
    Es = zeros(0)
    for line in eachline(joinpath("ref", @sprintf("S_%.1f__Jz_%.1f__Jxy_%.1f__G_%.1f__L_8.dat", 0.5S2, Jz, Jxy, G)))
        words = split(line)
        push!(Ts, parse(words[1]))
        push!(Es, parse(words[2]))
    end
    return Ts, Es
end

function E_QMC(T; S2=1, Jz=1.0, Jxy=1.0, G=0.0)
    p = Dict("Model"=>QuantumXXZ, "Lattice"=>chain_lattice,
             "S2"=>S2, "L"=>8,
             "Jz"=>Jz, "Jxy"=>Jxy,
             "G"=>G,
             "T"=>T,
             "MCS"=>8192,
            )
    res = runMC(p)
    return res["Energy"]
end

@testset "QuantumXXZ chain" begin
    const S2s = [1]
    const Js = [(1.0, 1.0), (1.0, 0.5), (1.0, 2.0),
                (0.0, 1.0),
                (-1.0, 1.0), (-1.0, 0.5), (-1.0, 2.0)]
    const Gs = [0.0, ]
    @testset "S=$(0.5S2)" for S2 in S2s
        @testset "Jz=$(J[1]), Jxy=$(J[2]) G=$G" for (G,J) in Iterators.product(Gs,Js)
            Ts, exacts = loaddata(S2=S2, Jz=J[1], Jxy=J[2], G=G)
            N = length(Ts)
            for (T,exact) in zip(Ts,exacts)
                ene = E_QMC(T, S2=S2, Jz=J[1], Jxy=J[2], G=G)
                if p_value(ene, exact) <= alpha/N
                    @show T
                    @show ene
                    @show exact
                end
                @test p_value(ene, exact) > alpha/N
            end
        end
    end
end
