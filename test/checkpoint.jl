@testset "checkpoint" begin
    @testset "$modelstr" for (modelstr, updatestrs) in (("Ising", ("local_update!", "SW_update!", "Wolff_update!")), 
                                                        ("Potts", ("local_update!", "SW_update!", "Wolff_update!")),
                                                        ("Clock", ("local_update!", "SW_update!", "Wolff_update!")),
                                                        ("XY", ("local_update!", "SW_update!", "Wolff_update!")),
                                                        ("QuantumXXZ", ("loop_update!",)),
                                                       )
        model = eval(Symbol(modelstr))
        @testset "$updatestr" for updatestr in updatestrs
            update = eval(Symbol(updatestr))
            p = Parameter("Model"=>model, "Lattice"=>"chain lattice", "L"=>8, "J"=>1.0, "T"=>1.0,
                     "Update Method" => update,
                     "Q"=>5, "S"=>0.5,
                     "Seed"=>SEED,
                     "MCS"=>100, "Thermalization"=>100,
                     "Checkpoint Interval"=>Inf)
            rm("cp_0.dat",force=true)
            runMC(p)
            p["MCS"] = 200
            p["Seed"] = 0
            res1 = runMC(p)
            rm("cp_0.dat",force=true)
            p["Seed"] = SEED
            res2 = runMC(p)
            p["Seed"] = 0
            res3 = runMC(p)
            for name in keys(res1)
                if name in ["MCS per Second", "Time per MCS"]
                    continue
                end
                @testset "$name" begin
                    @test mean(res1[name]) == mean(res2[name])
                    @test mean(res2[name]) == mean(res3[name])
                    @test stderror(res1[name]) == stderror(res2[name])
                    @test stderror(res2[name]) == stderror(res3[name])
                end
            end
            rm("cp_0.dat",force=true)
        end
    end
end
