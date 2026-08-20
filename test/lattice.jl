# Binds only the module name -- not `dim` -- so bare `dim` in this file keeps
# resolving to the one this package owns.
using Distributions: Distributions

@testset "lattice generation" begin
    @testset "interpolate does not mutate standard bravais definitions" begin
        lat = generatelattice(Parameter("Lattice" => "chain lattice", "L" => 4, "a" => 2.0))
        @test lat.latticevector == reshape([2.0], 1, 1)

        lat = generatelattice(Parameter("Lattice" => "chain lattice", "L" => 4, "a" => 3.0))
        @test lat.latticevector == reshape([3.0], 1, 1)

        generatelattice(Parameter("Lattice" => "triangular lattice", "L" => 2, "a" => 1.0))
        lat = generatelattice(Parameter("Lattice" => "honeycomb lattice", "L" => 2))
        @test lat.latticevector ≈ [sqrt(3.0) -sqrt(3.0)/2; 0.0 1.5]
    end

    @testset "dim is owned by this package" begin
        lat = generatelattice(Parameter("Lattice" => "square lattice", "L" => 4))
        @test SpinMonteCarlo.dim(lat) == 2
        @test SpinMonteCarlo.dim !== Distributions.dim
    end

    @testset "cubic lattice" begin
        lat = generatelattice(Parameter("Lattice" => "cubic lattice", "L" => [2, 2, 2]))
        @test dim(lat) == 3
        @test numsites(lat) == 8
        @test numbonds(lat) == 24
    end

    @testset "periodic boundary condition key" begin
        lat = generatelattice(Parameter("Lattice" => "square lattice",
                                        "L" => [3, 3],
                                        "Periodic Boundary Condition" => [false, false]))
        @test numbonds(lat) == 12

        @test_logs (:warn, r"Periodic Boudary Condition") begin
            lat = generatelattice(Parameter("Lattice" => "square lattice",
                                            "L" => [3, 3],
                                            "Periodic Boudary Condition" => [false, false]))
            @test numbonds(lat) == 12
        end

        lat = generatelattice(Parameter("Lattice" => "square lattice",
                                        "L" => [3, 3],
                                        "Periodic Boundary Condition" => [false, false],
                                        "Periodic Boudary Condition" => [true, true]))
        @test numbonds(lat) == 12
    end

    @testset "LatticeParameter indexing" begin
        lp = SpinMonteCarlo.LatticeParameter(SpinMonteCarlo.P(:dimension => 2))
        @test lp[:dimension] == 2

        lp[:foo] = "bar"
        @test lp[:foo] == "bar"
    end

    @testset "L parameter is not mutated" begin
        L = [4]
        generatelattice(Parameter("Lattice" => "square lattice", "L" => L))
        @test L == [4]
    end
end
