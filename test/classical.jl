const obsnames_ising = ["Energy", "Energy^2", "Specific Heat",
                        "|Magnetization|", "Magnetization^2", "Magnetization^4",
                        "Susceptibility", "Connected Susceptibility", "Binder Ratio"]
const obsnames_clock = ["Energy", "Energy^2", "Specific Heat",
                        "|Magnetization x|", "Magnetization x^2", "Magnetization x^4",
                        "Susceptibility x", "Connected Susceptibility x", "Binder Ratio x",
                        "|Magnetization y|", "Magnetization y^2", "Magnetization y^4",
                        "Susceptibility y", "Connected Susceptibility y", "Binder Ratio y",
                        "|Magnetization|", "|Magnetization|^2", "|Magnetization|^4",
                        "Susceptibility", "Connected Susceptibility", "Binder Ratio"]
function loaddata(filename, obsnames)
    Ts = zeros(0)
    res = Parameter(n => zeros(0) for n in obsnames)
    for line in eachline(filename)
        words = split(line)
        push!(Ts, parse(Float64, words[1]))
        for (i, n) in enumerate(obsnames)
            push!(res[n], parse(Float64, words[i + 1]))
        end
    end
    return Ts, res
end

function exact_ising_chain(L, J, T)
    Z = 0.0
    EZ = 0.0
    E2Z = 0.0
    for state in 0:(2 ^ L - 1)
        E = 0.0
        for i in 1:L
            si = (state >> (i - 1)) & 1 == 1 ? 1 : -1
            sj = (state >> (mod1(i + 1, L) - 1)) & 1 == 1 ? 1 : -1
            E -= J * si * sj
        end
        w = exp(-E / T)
        Z += w
        EZ += E * w
        E2Z += E^2 * w
    end
    e = EZ / Z / L
    e2 = E2Z / Z / L^2
    return Parameter("Energy" => e,
                     "Specific Heat" => L / T^2 * (e2 - e^2))
end

function exact_ashkin_teller_chain(L, Jsigma, Jtau, K, T)
    Z = 0.0
    EZ = 0.0
    for state in 0:(4 ^ L - 1)
        x = state
        sigma = zeros(Int, L)
        tau = zeros(Int, L)
        for i in 1:L
            sigma[i] = x & 1 == 1 ? 1 : -1
            x >>= 1
            tau[i] = x & 1 == 1 ? 1 : -1
            x >>= 1
        end

        E = 0.0
        for i in 1:L
            j = mod1(i + 1, L)
            E -= Jsigma * sigma[i] * sigma[j]
            E -= Jtau * tau[i] * tau[j]
            E -= K * sigma[i] * sigma[j] * tau[i] * tau[j]
        end
        w = exp(-E / T)
        Z += w
        EZ += E * w
    end
    return EZ / Z / L
end

@testset "chain transfer matrix regressions" begin
    Ts = [0.5, 1.0, 2.0]
    L = 8
    nT = length(Ts)

    @testset "Ising chain J=$J $upstr" for (J, updates) in
                                           [(1.0,
                                             ("local_update!",
                                              "SW_update!",
                                              "Wolff_update!")),
                                            (-1.0,
                                             ("local_update!",
                                              "SW_update!",
                                              "Wolff_update!"))],
                                           upstr in updates

        p = Parameter("Model" => Ising,
                      "Lattice" => "chain lattice",
                      "L" => L,
                      "J" => J,
                      "MCS" => MCS,
                      "Thermalization" => Therm,
                      "Update Method" => eval(Symbol(upstr)))
        res1 = []
        res2 = []
        for T in Ts
            p["T"] = T
            p["Seed"] = SEED
            push!(res1, runMC(p))
            p["Seed"] = SEED2
            push!(res2, runMC(p))
        end

        @testset "$obsname" for obsname in ("Energy", "Specific Heat")
            for i in 1:nT
                exact = exact_ising_chain(L, J, Ts[i])[obsname]
                mc1 = res1[i][obsname]
                mc2 = res2[i][obsname]
                if !(p_value(mc1, exact) > alpha / nT ||
                     p_value(mc2, exact) > alpha / nT)
                    @show J
                    @show upstr
                    @show Ts[i]
                    @show obsname
                    @show exact
                    @show mc1, p_value(mc1, exact)
                    @show mc2, p_value(mc2, exact)
                end
                @test p_value(mc1, exact) > alpha / nT ||
                      p_value(mc2, exact) > alpha / nT
            end
        end
    end

    @testset "AshkinTeller chain local energy" begin
        L = 4
        Jsigma = 1.0
        Jtau = 0.8
        K = 0.3
        T = 1.5
        exact = exact_ashkin_teller_chain(L, Jsigma, Jtau, K, T)
        p = Parameter("Model" => AshkinTeller,
                      "Lattice" => "chain lattice",
                      "L" => L,
                      "Jsigma" => Jsigma,
                      "Jtau" => Jtau,
                      "K" => K,
                      "T" => T,
                      "MCS" => MCS,
                      "Thermalization" => Therm,
                      "Update Method" => local_update!)
        p["Seed"] = SEED
        res1 = runMC(p)
        p["Seed"] = SEED2
        res2 = runMC(p)
        @test p_value(res1["Energy"], exact) > alpha ||
              p_value(res2["Energy"], exact) > alpha
    end
end

@testset "cluster update regressions" begin
    @testset "Ising SWInfo clustermag" begin
        lat = generatelattice(Parameter("Lattice" => "chain lattice", "L" => 8))
        model = Ising(lat, SEED)
        preflip_magnetization = sum(model.spins)
        sw = SW_update!(model, 1.0, fill(-1.0, SpinMonteCarlo.numbondtypes(model)))
        @test hasproperty(sw, :clustermag)
        @test sum(sw.clustermag) == preflip_magnetization
        @test length(sw.clustermag) == SpinMonteCarlo.numclusters(sw)
    end

    @testset "Potts AF cluster guard" begin
        lat = generatelattice(Parameter("Lattice" => "chain lattice", "L" => 4))
        model = Potts(lat, 3, SEED)
        @test_throws ArgumentError SW_update!(model, 1.0, [-1.0])
        @test_throws ArgumentError Wolff_update!(model, 1.0, [-1.0])
    end
end

@testset "helicity modulus regressions" begin
    @testset "$modelstr J scaling" for modelstr in ("XY", "Clock")
        lat = generatelattice(Parameter("Lattice" => "square lattice", "L" => 4))
        model = modelstr == "XY" ? XY(lat, SEED) : Clock(lat, 6, SEED)
        model.spins .= 1
        r1 = simple_estimator(model, 1.0, fill(1.0, SpinMonteCarlo.numbondtypes(model)))
        r2 = simple_estimator(model, 1.0, fill(2.0, SpinMonteCarlo.numbondtypes(model)))
        @test r2["Helicity Modulus x"] / r1["Helicity Modulus x"] ≈ 2.0
        @test r2["Helicity Modulus y"] / r1["Helicity Modulus y"] ≈ 2.0
    end

    @testset "XY cubic z helicity" begin
        p = Parameter("Model" => XY,
                      "Lattice" => "cubic lattice",
                      "L" => [2, 2, 2],
                      "J" => 1.0,
                      "T" => 1.0,
                      "MCS" => 16,
                      "Thermalization" => 16,
                      "Seed" => SEED,
                      "Update Method" => local_update!)
        res = runMC(p)
        @test haskey(res, "Helicity Modulus z")
    end
end

function parse_filename(filename, ::Union{Type{Ising},Type{XY}})
    m = match(r"^J_([\d.-]*)__N_([\d.-]*).dat$", filename)
    if m === nothing
        return nothing
    end
    p = Parameter()
    p["J"] = parse(Float64, m.captures[1])
    p["N"] = parse(Int, m.captures[2])
    return p
end

function parse_filename(filename, ::Union{Type{Potts},Type{Clock}})
    m = match(r"^Q_(\d*)__J_([\d.-]*)__N_(\d*).dat$", filename)
    if m === nothing
        return nothing
    end
    p = Parameter()
    p["Q"] = parse(Int, m.captures[1])
    p["J"] = parse(Float64, m.captures[2])
    p["N"] = parse(Int, m.captures[3])
    return p
end

@testset "$modelstr" for (modelstr, pnames, obsnames) in
                         [("Ising", ("J", "N"), obsnames_ising),
                          ("Potts", ("Q", "J", "N"), obsnames_ising),
                          ("Clock", ("Q", "J", "N"), obsnames_clock),
                          ("XY", ("J", "N"), obsnames_clock)]
    model = eval(Symbol(modelstr))
    for filename in readdir(joinpath("ref", modelstr))
        p = parse_filename(filename, model)
        if p === nothing
            continue
        end
        testname = ""
        for pname in pnames
            testname *= "$(pname)=$(p[pname]) "
        end
        @testset "$testname" begin
            Ts, exacts = loaddata(joinpath("ref", modelstr, filename), obsnames)
            nT = length(Ts)
            p["Model"] = model
            p["Lattice"] = "fully connected graph"
            p["J"] = p["J"] / p["N"]
            p["MCS"] = MCS
            p["Thermalization"] = Therm
            @testset "$(upstr)" for upstr in ("local_update!",
                                              "SW_update!",
                                              "Wolff_update!")
                p["Update Method"] = eval(Symbol(upstr))
                res1 = []
                res2 = []
                for i in 1:nT
                    p["T"] = Ts[i]
                    p["Seed"] = SEED
                    push!(res1, runMC(p))
                    p["Seed"] = SEED2
                    push!(res2, runMC(p))
                end
                @testset "$n" for n in obsnames
                    for i in 1:nT
                        T = Ts[i]
                        exact = exacts[n][i]
                        r1 = res1[i]
                        r2 = res2[i]
                        ## single MC test may fail.
                        mc1 = r1[n]
                        mc2 = r2[n]
                        ex = exact
                        if !(p_value(mc1, exact) > alpha / nT ||
                             p_value(mc2, exact) > alpha / nT)
                            @show T
                            @show exact
                            @show mc1, p_value(mc1, exact)
                            @show mc2, p_value(mc2, exact)
                        end
                        @test p_value(mc1, exact) > alpha / nT ||
                              p_value(mc2, exact) > alpha / nT
                    end
                end
            end
        end
    end
end
