using Random

@testset "RNG construction" begin
    lat = generatelattice(Parameter("Lattice" => "chain lattice", "L" => 4))

    @testset "model constructors consume deterministic spin initialization RNG" begin
        ising_rng = Xoshiro(SEED)
        @test Ising(lat, SEED).spins == rand(ising_rng, [1, -1], 1, numsites(lat))

        Q_potts = 3
        potts_rng = Xoshiro(SEED)
        @test Potts(lat, Q_potts, SEED).spins ==
              rand(potts_rng, 1:Q_potts, 1, numsites(lat))

        Q_clock = 4
        clock_rng = Xoshiro(SEED)
        @test Clock(lat, Q_clock, SEED).spins ==
              rand(clock_rng, 1:Q_clock, 1, numsites(lat))

        xy_rng = Xoshiro(SEED)
        @test XY(lat, SEED).spins == rand(xy_rng, 1, numsites(lat))

        ashkin_teller_rng = Xoshiro(SEED)
        @test AshkinTeller(lat, SEED).spins ==
              rand(ashkin_teller_rng, [1, -1], 2, numsites(lat))

        S = 1 // 2
        S2 = round(Int, 2S)
        quantum_xxz_rng = Xoshiro(SEED)
        @test QuantumXXZ(lat, S, SEED).spins ==
              rand(quantum_xxz_rng, [1, -1], 1, numsites(lat) * S2)
    end

    @testset "runMC(param) preserves constructor RNG state" begin
        param = Parameter("Model" => Ising,
                          "Lattice" => "chain lattice",
                          "L" => 4,
                          "J" => 1.0,
                          "T" => 1.5,
                          "MCS" => 200,
                          "Thermalization" => 100,
                          "Seed" => SEED,
                          "Update Method" => local_update!)
        direct = runMC(param)
        model = param["Model"](param)
        explicit = runMC(model, param)

        @test mean(direct["Energy"]) == mean(explicit["Energy"])
        @test mean(direct["|Magnetization|"]) == mean(explicit["|Magnetization|"])
        @test mean(direct["Magnetization^2"]) == mean(explicit["Magnetization^2"])
    end

    @testset "param RNG type switches model type and stream" begin
        default_param = Parameter("Model" => Ising,
                                  "Lattice" => "chain lattice",
                                  "L" => 4,
                                  "J" => 1.0,
                                  "T" => 1.5,
                                  "MCS" => 200,
                                  "Thermalization" => 100,
                                  "Seed" => SEED,
                                  "Update Method" => local_update!)
        mt_param = copy(default_param)
        mt_param["RNG"] = MersenneTwister

        @test typeof(default_param["Model"](default_param)) == Ising{Xoshiro}
        @test typeof(mt_param["Model"](mt_param)) == Ising{MersenneTwister}

        default_result = runMC(default_param)
        mt_result = runMC(mt_param)
        @test mean(default_result["Energy"]) != mean(mt_result["Energy"]) ||
              mean(default_result["|Magnetization|"]) != mean(mt_result["|Magnetization|"])
    end
end
