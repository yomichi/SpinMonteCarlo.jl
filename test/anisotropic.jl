@testset "anisotropic" begin
    param = Dict{String,Any}("Model"=>Ising, "UpdateMethod"=>Wolff_update!,
                             "T"=>1.0, "L"=>8,
                             "MCS"=>8192, "Thermalization"=>8192,
                            )
    @testset "dimer-chain" begin
        srand(SEED)
        param["Lattice"] = dimer_lattice
        param["J"] = [1.0]
        ref = runMC(param)

        param["Lattice"] = chain_lattice
        param["J"] = [1.0, 0.0]
        obs = runMC(param)
        diff = obs["Energy"] - ref["Energy"]
        @test mean(diff) < confidence_interval(diff, 0.99)
        param["J"] = [0.0, 1.0]
        obs = runMC(param)
        diff = obs["Energy"] - ref["Energy"]
        @test mean(diff) < confidence_interval(diff, 0.99)
    end
    @testset "chain-square" begin
        srand(SEED)
        param["Lattice"] = chain_lattice
        param["J"] = [1.0, 1.0]
        ref = runMC(param)
        param["Lattice"] = square_lattice
        param["J"] = [1.0, 0.0]
        obs = runMC(param)
        diff = obs["Energy"] - ref["Energy"]
        @test mean(diff) < confidence_interval(diff, 0.99)

        param["J"] = [0.0, 1.0]
        obs = runMC(param)
        diff = obs["Energy"] - ref["Energy"]
        @test mean(diff) < confidence_interval(diff, 0.99)
    end

    @testset "chain-triangular" begin
        srand(SEED)
        param["Lattice"] = chain_lattice
        param["J"] = [1.0, 1.0]
        ref = runMC(param)
        param["Lattice"] = triangular_lattice
        param["J"] = [1.0, 0.0, 0.0]
        obs = runMC(param)
        diff = obs["Energy"] - ref["Energy"]
        @test mean(diff) < confidence_interval(diff, 0.99)
        param["J"] = [0.0, 1.0, 0.0]
        obs = runMC(param)
        diff = obs["Energy"] - ref["Energy"]
        @test mean(diff) < confidence_interval(diff, 0.99)
        param["J"] = [0.0, 0.0, 1.0]
        obs = runMC(param)
        diff = obs["Energy"] - ref["Energy"]
        @test mean(diff) < confidence_interval(diff, 0.99)
    end

    @testset "chain-cubic" begin
        srand(SEED)
        param["Lattice"] = chain_lattice
        param["J"] = [1.0, 1.0]
        ref = runMC(param)
        param["Lattice"] = cubic_lattice

        param["J"] = [1.0, 0.0, 0.0]
        obs = runMC(param)
        diff = obs["Energy"] - ref["Energy"]
        @test mean(diff) < confidence_interval(diff, 0.99)

        param["J"] = [0.0, 1.0, 0.0]
        obs = runMC(param)
        diff = obs["Energy"] - ref["Energy"]
        @test mean(diff) < confidence_interval(diff, 0.99)

        param["J"] = [0.0, 0.0, 1.0]
        obs = runMC(param)
        diff = obs["Energy"] - ref["Energy"]
        @test mean(diff) < confidence_interval(diff, 0.99)
    end

    @testset "square-cubic" begin
        srand(SEED)
        param["Lattice"] = square_lattice
        param["J"] = [1.0, 1.0]
        ref = runMC(param)
        param["Lattice"] = cubic_lattice

        param["J"] = [1.0, 1.0, 0.0]
        obs = runMC(param)
        diff = obs["Energy"] - ref["Energy"]
        @test mean(diff) < confidence_interval(diff, 0.99)
        param["J"] = [1.0, 0.0, 1.0]
        obs = runMC(param)
        diff = obs["Energy"] - ref["Energy"]
        @test mean(diff) < confidence_interval(diff, 0.99)
        param["J"] = [0.0, 1.0, 1.0]
        obs = runMC(param)
        diff = obs["Energy"] - ref["Energy"]
        @test mean(diff) < confidence_interval(diff, 0.99)
    end
end
