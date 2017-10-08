using SpecialFunctions

function energy_ising_dimer(T)
    beta = 1.0/T
    Z = exp(beta) + exp(-beta)
    return 0.5(exp(-beta) - exp(beta))/Z
end

function energy_potts_dimer(Q,T)
    eb = exp(1.0/T)
    return -0.5eb/(eb+Q-1)
end

function energy_xy_dimer(T)
    return -0.5besseli(1,1.0/T) / besseli(0, 1.0/T)
end

function energy_clock_dimer(Q,T)
    beta = 1.0/T
    Z = 0.0
    E = 0.0
    for q in 1:Q
        ene = -cospi(2*(q-1)/Q)
        z = exp(-ene*beta)
        Z += z
        E += ene*z
    end
    return 0.5E/Z
end

function energy_dimer(param::Dict)
    model = param["Model"]
    if model == Ising
        return energy_ising_dimer(param["T"])
    elseif model == Potts
        return energy_potts_dimer(param["Q"], param["T"])
    elseif model == XY
        return energy_xy_dimer(param["T"])
    elseif model == Clock
        return energy_clock_dimer(param["Q"], param["T"])
    else
        error("unknown model")
    end
end

@testset "dimer energy" begin
    param = Dict{String,Any}("Model" => Ising, "Lattice" => dimer_lattice,
                              "J" => [1.0],
                              "MCS" => 32768, "Thermalization" => 0,
                             )
    const Ts = [0.3, 1.0, 3.0]
    @testset "$upname" for (upname, method) in [("local update", local_update!),
                                                ("Swendsen-Wang", SW_update!),
                                                ("Wolff", Wolff_update!),
                                               ]
        param["UpdateMethod"] = method
        @testset "$modelname" for (modelname, model) in [("Ising", Ising),
                                                         ("Potts", Potts),
                                                         ("XY", XY),
                                                         ("Clock", Clock),
                                                        ]
            srand(SEED)
            param["Model"] = model
            if model == Potts
                param["Q"] = 3
            elseif model == Clock
                param["Q"] = 5
            end
            @testset "T = $T" for T in Ts
                param["T"] = T
                obs = runMC(param)
                exact = energy_dimer(param)
                @test abs(mean(obs["Energy"]) - exact) < 3.0*stderror(obs["Energy"])
            end
        end
    end
end
