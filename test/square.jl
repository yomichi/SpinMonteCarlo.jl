@testset "square lattice" begin

    @testset "$mname" for (mname, model, obsnames) in [("Ising", Ising, obsnames_ising),
                                                       ("Potts", Potts, obsnames_potts),
                                                       ("XY", XY, obsnames_xy),
                                                       ("Clock", Clock, obsnames_clock),
                                                      ]
        srand(SEED)
        param = Dict{String,Any}("Model" => model, "Lattice" => square_lattice, "J" => 1.0,
                                  "L" => 6, "T" => 3.0,
                                  "UpdateMethod" => local_update!,
                                  "MCS" => MCS, "Thermalization" => Therm,
                                 )
        if model==Potts
            param["Q"] = 3
        elseif model==Clock
            param["Q"] = 5
        end
        localupdate = runMC(param)
        param["UpdateMethod"] = SW_update!
        sw = runMC(param)
        param["UpdateMethod"] = Wolff_update!
        wolff = runMC(param)

        @testset "$name" for name in obsnames
            diff = localupdate[name] - sw[name]
            @test abs(mean(diff)) < confidence_interval(diff, conf_ratio)
            diff = localupdate[name] - wolff[name]
            @test abs(mean(diff)) < confidence_interval(diff, conf_ratio)
        end
    end
end
