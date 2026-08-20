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

    @testset "childseed fixed values" begin
        @test SpinMonteCarlo.childseed(0, 0) == 5197578548964807871
        @test SpinMonteCarlo.childseed(137, 1) == 0x8e7e968bcd197b0e
        @test SpinMonteCarlo.childseed(137, 2) == 0x1c23f3da2a9bedd6
        @test SpinMonteCarlo.childseed(137, 3) == 0xaa13ff27c06fbe76
        @test SpinMonteCarlo.childseed(-2, 5) == 2745338216157443808
        @test SpinMonteCarlo.childseed(0, -5) == 741713033729228089
        @test SpinMonteCarlo.childseed(-2, 5) != SpinMonteCarlo.childseed(0, -5)
        @test SpinMonteCarlo.childseed(0, 0) != 0
    end

    @testset "makerng derives child streams from IDs" begin
        rngsample(rng) = rand(rng, UInt64, 4)

        seeded = Parameter("Seed" => SEED)
        @test rngsample(SpinMonteCarlo.makerng(seeded)) == rngsample(Xoshiro(SEED))

        with_id = Parameter("Seed" => SEED, "ID" => 1)
        @test rngsample(SpinMonteCarlo.makerng(with_id)) != rngsample(Xoshiro(SEED))
        @test rngsample(SpinMonteCarlo.makerng(with_id)) ==
              rngsample(SpinMonteCarlo.makerng(Parameter("Seed" => SEED, "ID" => 1)))

        other_id = Parameter("Seed" => SEED, "ID" => 2)
        @test rngsample(SpinMonteCarlo.makerng(with_id)) !=
              rngsample(SpinMonteCarlo.makerng(other_id))

        mt_param = Parameter("Seed" => SEED, "ID" => 1, "RNG" => MersenneTwister)
        @test SpinMonteCarlo.makerng(mt_param) isa MersenneTwister
    end

    @testset "runMC array derives per-ID child streams" begin
        function make_params()
            return [Parameter("Model" => Ising,
                              "Lattice" => "chain lattice",
                              "L" => 4,
                              "J" => 1.0,
                              "T" => 1.5,
                              "MCS" => 200,
                              "Thermalization" => 100,
                              "Seed" => SEED,
                              "Update Method" => local_update!) for _ in 1:3]
        end
        result_signature(res) = (mean(res["Energy"]),
                                 mean(res["|Magnetization|"]),
                                 mean(res["Magnetization^2"]))

        signatures1 = result_signature.(runMC(make_params()))
        signatures2 = result_signature.(runMC(make_params()))
        @test length(unique(signatures1)) == length(signatures1)
        @test signatures1 == signatures2
    end

    @testset "parallel runMC array smoke test" begin
        params = [Parameter("Model" => Ising,
                            "Lattice" => "chain lattice",
                            "L" => 4,
                            "J" => 1.0,
                            "T" => 1.5,
                            "MCS" => 100,
                            "Thermalization" => 50,
                            "Seed" => SEED,
                            "Update Method" => local_update!) for _ in 1:2]
        results = runMC(params; parallel=true)
        @test length(results) == 2
        @test all(haskey(result, "Energy") for result in results)
    end

    @testset "child seed derivation requires an integer seed" begin
        seeded = Parameter("Lattice" => "chain lattice", "L" => 4,
                           "Seed" => "not an integer", "ID" => 1)
        @test_throws ArgumentError Ising(seeded)

        # Without an "ID" nothing is derived, so the seed goes straight to the RNG
        # constructor and whatever Julia accepts there keeps working.
        undivided = Parameter("Lattice" => "chain lattice", "L" => 4,
                              "Seed" => "not an integer")
        @test Ising(undivided) isa Ising{Xoshiro}
    end

    @testset "RNG types need only the constructor the parameters call for" begin
        # makerng calls T(seed) when "Seed" is given and T() when it is not, so a
        # generator that offers only T() is usable as long as no seed is supplied.
        for T in (Random.RandomDevice, Random.TaskLocalRNG)
            seedless = Parameter("Lattice" => "chain lattice", "L" => 4, "RNG" => T)
            @test Ising(seedless) isa Ising{T}

            seeded = copy(seedless)
            seeded["Seed"] = SEED
            @test_throws MethodError Ising(seeded)
        end
    end
end
