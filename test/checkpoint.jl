@testset "checkpoint" begin
    @testset "$modelstr" for modelstr in ("Ising",
                                          #"Potts",
                                          "Clock",
                                          #"XY",
                                          #"QuantumXXZ",
                                         )
        model = eval(Symbol(modelstr))
        p = Dict("Model"=>model, "Lattice"=>chain_lattice, "L"=>8, "J"=>1.0, "T"=>1.0,
                 "Q"=>5, "S"=>0.5,
                 "Checkpoint Interval"=>Inf)
        rm("cp_1.jld2",force=true)
        res1 = runMC(p)
        res2 = runMC(p)
        @testset "$name" for name in keys(res1)
            @test mean(res1[name]) == mean(res2[name])
            @test stderror(res1[name]) == stderror(res2[name])
        end
        rm("cp_1.jld2",force=true)
    end
end
