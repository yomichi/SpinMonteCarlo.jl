@testset "square lattice" begin

    @testset "Ising" begin
        srand(SEED)
        param = Dict{String,Any}("Model" => Ising, "Lattice" => square_lattice, "J" => [1.0,1.0],
                                  "L" => 6, "T" => 3.0,
                                  "UpdateMethod" => local_update!,
                                  "MCS" => 8192, "Thermalization" => 8192,
                                 )
        localupdate = runMC(param)
        param["UpdateMethod"] = SW_update!
        sw = runMC(param)
        param["UpdateMethod"] = Wolff_update!
        wolff = runMC(param)

        @testset "$name" for name in obsnames_ising
            diff = localupdate[name] - sw[name]
            @test abs(mean(diff)) < confidence_interval(diff, conf_ratio)
            diff = localupdate[name] - wolff[name]
            @test abs(mean(diff)) < confidence_interval(diff, conf_ratio)
        end
    end

    @testset "Potts" begin
        srand(SEED)
        param = Dict{String,Any}("Model" => Potts, "Lattice" => square_lattice, "J" => [1.0,1.0],
                                  "Q" => 3, "L" => 6, "T" => 2.0,
                                  "UpdateMethod" => local_update!,
                                  "MCS" => 8192, "Thermalization" => 8192,
                                 )
        localupdate = runMC(param)
        param["UpdateMethod"] = SW_update!
        sw = runMC(param)
        param["UpdateMethod"] = Wolff_update!
        wolff = runMC(param)

        @testset "$name" for name in obsnames_potts
            diff = localupdate[name] - sw[name]
            @test abs(mean(diff)) < confidence_interval(diff, conf_ratio)
            diff = localupdate[name] - wolff[name]
            @test abs(mean(diff)) < confidence_interval(diff, conf_ratio)
        end
    end

    @testset "XY" begin
        srand(SEED)
        param = Dict{String,Any}("Model" => XY, "Lattice" => square_lattice, "J" => [1.0,1.0],
                                  "L" => 6, "T" => 10.0,
                                  "UpdateMethod" => local_update!,
                                  "MCS" => 8192, "Thermalization" => 8192,
                                 )
        localupdate = runMC(param)
        param["UpdateMethod"] = SW_update!
        sw = runMC(param)
        param["UpdateMethod"] = Wolff_update!
        wolff = runMC(param)

        @testset "$name" for name in obsnames_xy
            diff = localupdate[name] - sw[name]
            @test abs(mean(diff)) < confidence_interval(diff, conf_ratio)
            diff = localupdate[name] - wolff[name]
            @test abs(mean(diff)) < confidence_interval(diff, conf_ratio)
        end
    end

    @testset "Clock" begin
        srand(SEED)
        param = Dict{String,Any}("Model" => Clock, "Lattice" => square_lattice, "J" => [1.0,1.0],
                                  "Q" => 5, "L" => 6, "T" => 5.0,
                                  "UpdateMethod" => local_update!,
                                  "MCS" => 8192, "Thermalization" => 8192,
                                 )
        localupdate = runMC(param)
        param["UpdateMethod"] = SW_update!
        sw = runMC(param)
        param["UpdateMethod"] = Wolff_update!
        wolff = runMC(param)

        @testset "$name" for name in obsnames_xy
            diff = localupdate[name] - sw[name]
            @test abs(mean(diff)) < confidence_interval(diff, conf_ratio)
            diff = localupdate[name] - wolff[name]
            @test abs(mean(diff)) < confidence_interval(diff, conf_ratio)
        end
    end
end
