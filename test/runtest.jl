using Base.Test

include("../src/SpinMonteCarlo.jl")
using SpinMonteCarlo

const obsnames_ising = ["Magnetization", "|Magnetization|", "Magnetization^2", "Magnetization^4", "Binder Ratio",
                        "Susceptibility", "Connected Susceptibility",
                        "Energy", "Energy^2", "Specific Heat",
                       ]
const obsnames_potts = obsnames_ising

const obsnames_xy = ["|Magnetization|", "|Magnetization|^2", "|Magnetization|^4",
                        "Binder Ratio", "Susceptibility", "Connected Susceptibility",
                        "Magnetization x", "|Magnetization x|", "Magnetization x^2", "Magnetization x^4",
                        "Binder Ratio x", "Susceptibility x", "Connected Susceptibility x",
                        "Magnetization y", "|Magnetization y|", "Magnetization y^2", "Magnetization y^4",
                        "Binder Ratio y", "Susceptibility y", "Connected Susceptibility y",
                        "Helicity Modulus x", "Helicity Modulus y",
                        "Energy", "Energy^2", "Specific Heat",
                       ]
const obsnames_clock = obsnames_xy

@testset "Ising" begin
    srand(1234567)
    param = Dict{String,Any}("Model" => Ising, "Lattice" => square_lattice,
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
        @test abs(mean(localupdate[name]) - mean(sw[name])) < 3.0(stderror(localupdate[name])+stderror(sw[name]))
        @test abs(mean(localupdate[name]) - mean(wolff[name])) < 3.0(stderror(localupdate[name])+stderror(wolff[name]))
    end
end

@testset "Potts" begin
    srand(1234567)
    param = Dict{String,Any}("Model" => Potts, "Lattice" => square_lattice,
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
        @test abs(mean(localupdate[name]) - mean(sw[name])) < 3.0(stderror(localupdate[name])+stderror(sw[name]))
        @test abs(mean(localupdate[name]) - mean(wolff[name])) < 3.0(stderror(localupdate[name])+stderror(wolff[name]))
    end
end

@testset "XY" begin
    srand(1234567)
    param = Dict{String,Any}("Model" => XY, "Lattice" => square_lattice,
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
        @test abs(mean(localupdate[name]) - mean(sw[name])) < 3.0(stderror(localupdate[name])+stderror(sw[name]))
        @test abs(mean(localupdate[name]) - mean(wolff[name])) < 3.0(stderror(localupdate[name])+stderror(wolff[name]))
    end
end

@testset "Clock" begin
    srand(1234567)
    param = Dict{String,Any}("Model" => Clock, "Lattice" => square_lattice,
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
        @test abs(mean(localupdate[name]) - mean(sw[name])) < 3.0(stderror(localupdate[name])+stderror(sw[name]))
        @test abs(mean(localupdate[name]) - mean(wolff[name])) < 3.0(stderror(localupdate[name])+stderror(wolff[name]))
    end
end
