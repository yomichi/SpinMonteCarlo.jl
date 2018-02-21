using SpecialFunctions

macro dimer_test(o, name, exactval, conf_ratio)
    quote
        c = confidence_interval($(esc(o))[$(esc(name))], $(esc(conf_ratio)))
        if c == 0.0
            c = sqrt(nextfloat($(esc(exactval))) - $(esc(exactval)))
        end
        @test abs(mean($(esc(o))[$(esc(name))]) - $(esc(exactval))) <= c
    end
end

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

function energy_TFI_dimer(T,J,G)
    sz = [0.5 0.0; 0.0 -0.5]
    sx = [0.0 0.5; 0.5 0.0]
    H = J.*kron(sz,sz) .- G.*(kron(sx,eye(2)) .+ kron(eye(2),sx))
    rho = expm((-1.0/T).*H)
    return trace(H*rho)/trace(rho)
end

function mag2_TFI_dimer(T,J,G)
    sz = [0.5 0.0; 0.0 -0.5]
    sx = [0.0 0.5; 0.5 0.0]
    H = J.*kron(sz,sz) .- G.*(kron(sx,eye(2)) .+ kron(eye(2),sx))
    rho = expm((-1.0/T).*H)
    sz_tot = 0.5.*(kron(sz, eye(2)) + kron(eye(2), sz))
    return trace(sz_tot*sz_tot*rho)/trace(rho)
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
    elseif model == TransverseFieldIsing
        return energy_TFI_dimer(param["T"], param["J"], param["Gamma"])
    else
        error("unknown model")
    end
end

@testset "dimer energy" begin
    param = Dict{String,Any}("Model" => Ising, "Lattice" => dimer_lattice,
                              "J" => 1.0,
                              "MCS" => 10000, "Thermalization" => 1000,
                             )
    const Ts = [0.1, 0.5, 1.0]
    @testset "$upname" for (upname, method) in [
                                                ("local update", local_update!),
                                                ("Swendsen-Wang", SW_update!),
                                                ("Wolff", Wolff_update!),
                                               ]
        param["UpdateMethod"] = method
        @testset "$modelname" for (modelname, model) in [
                                                         ("Ising", Ising),
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
                @test abs(mean(obs["Energy"]) - exact) <= confidence_interval(obs["Energy"],conf_ratio)
            end
        end
    end
    @testset "TransverseFieldIsing" begin
        srand(SEED)
        Js = [1.0]
        Gs = [0.0, 0.1, 0.5, 1.0, 5.0]
        param["Model"] = TransverseFieldIsing
        @testset "J = $J, G = $G, T = $T" for (J,G,T) in Iterators.product(Js,Gs,Ts)
            param["J"] = J
            param["Gamma"] = G
            param["T"] = T
            obs = runMC(param)
            exact_energy = energy_dimer(param)
            exact_mag2 = mag2_TFI_dimer(T,J,G)
            @dimer_test(obs, "Energy", exact_energy, conf_ratio)
            @dimer_test(obs, "Magnetization", 0.0, conf_ratio)
            @dimer_test(obs, "Magnetization^2", exact_mag2, conf_ratio)
        end
    end
end
