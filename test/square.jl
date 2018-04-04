@testset "square lattice" begin

    @testset "$mname" for (mname, model, obsnames) in [("Ising", Ising, obsnames_ising),
                                                       ("Potts", Potts, obsnames_potts),
                                                       ("XY", XY, obsnames_xy),
                                                       ("Clock", Clock, obsnames_clock),
                                                      ]
        srand(SEED)
        param = Dict{String,Any}("Model" => model, "Lattice" => square_lattice, "J" => 1.0,
                                  "L" => 6, "T" => 3.5,
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
            if p_value(localupdate[name], sw[name]) <= alpha
                @show localupdate[name], sw[name]
            end
            @test p_value(localupdate[name], sw[name]) > alpha

            if p_value(localupdate[name], wolff[name]) <= alpha
                @show localupdate[name], wolff[name]
            end
            @test p_value(localupdate[name], wolff[name]) > alpha

            if p_value(sw[name], wolff[name]) <= alpha
                @show sw[name], wolff[name]
            end
            @test p_value(sw[name], wolff[name]) > alpha
        end
    end
end
