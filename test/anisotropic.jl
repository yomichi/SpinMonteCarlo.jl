@testset "anisotropic" begin
    param = Dict{String,Any}("Model"=>Ising, "UpdateMethod"=>Wolff_update!,
                             "T"=>1.0, "L"=>8,
                             "MCS"=>MCS, "Thermalization"=>Therm,
                            )
    @testset "dimer-chain" begin
        srand(SEED)
        param["Lattice"] = dimer_lattice
        param["J"] = [1.0]
        ref = runMC(param)

        param["Lattice"] = chain_lattice
        param["J"] = [1.0, 0.0]
        obs = runMC(param)
        @test p_value(obs["Energy"], ref["Energy"]) > alpha
        param["J"] = [0.0, 1.0]
        obs = runMC(param)
        @test p_value(obs["Energy"], ref["Energy"]) > alpha
    end
    @testset "chain-square" begin
        srand(SEED)
        param["Lattice"] = chain_lattice
        param["J"] = [1.0, 1.0]
        ref = runMC(param)
        param["Lattice"] = square_lattice
        param["J"] = [1.0, 0.0]
        obs = runMC(param)
        @test p_value(obs["Energy"], ref["Energy"]) > alpha

        param["J"] = [0.0, 1.0]
        obs = runMC(param)
        @test p_value(obs["Energy"], ref["Energy"]) > alpha
    end

    @testset "chain-triangular" begin
        srand(SEED)
        param["Lattice"] = chain_lattice
        param["J"] = [1.0, 1.0]
        ref = runMC(param)
        param["Lattice"] = triangular_lattice
        param["J"] = [1.0, 0.0, 0.0]
        obs = runMC(param)
        @test p_value(obs["Energy"], ref["Energy"]) > alpha
        param["J"] = [0.0, 1.0, 0.0]
        obs = runMC(param)
        @test p_value(obs["Energy"], ref["Energy"]) > alpha
        param["J"] = [0.0, 0.0, 1.0]
        obs = runMC(param)
        @test p_value(obs["Energy"], ref["Energy"]) > alpha
    end

    @testset "chain-cubic" begin
        srand(SEED)
        param["Lattice"] = chain_lattice
        param["J"] = [1.0, 1.0]
        ref = runMC(param)
        param["Lattice"] = cubic_lattice

        param["J"] = [1.0, 0.0, 0.0]
        obs = runMC(param)
        @test p_value(obs["Energy"], ref["Energy"]) > alpha

        param["J"] = [0.0, 1.0, 0.0]
        obs = runMC(param)
        @test p_value(obs["Energy"], ref["Energy"]) > alpha

        param["J"] = [0.0, 0.0, 1.0]
        obs = runMC(param)
        @test p_value(obs["Energy"], ref["Energy"]) > alpha
    end

    @testset "square-cubic" begin
        srand(SEED)
        param["Lattice"] = square_lattice
        param["J"] = [1.0, 1.0]
        ref = runMC(param)
        param["Lattice"] = cubic_lattice

        param["J"] = [1.0, 1.0, 0.0]
        obs = runMC(param)
        @test p_value(obs["Energy"], ref["Energy"]) > alpha
        param["J"] = [1.0, 0.0, 1.0]
        obs = runMC(param)
        @test p_value(obs["Energy"], ref["Energy"]) > alpha
        param["J"] = [0.0, 1.0, 1.0]
        obs = runMC(param)
        @test p_value(obs["Energy"], ref["Energy"]) > alpha
    end
end
