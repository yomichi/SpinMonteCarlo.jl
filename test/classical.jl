const obsnames_ising = ["Energy", "Energy^2", "Specific Heat",
                        "|Magnetization|", "Magnetization^2", "Magnetization^4",
                        "Susceptibility", "Connected Susceptibility", "Binder Ratio",
                       ]
const obsnames_clock = ["Energy", "Energy^2", "Specific Heat",
                        "|Magnetization x|", "Magnetization x^2", "Magnetization x^4",
                        "Susceptibility x", "Connected Susceptibility x", "Binder Ratio x",
                        "|Magnetization y|", "Magnetization y^2", "Magnetization y^4",
                        "Susceptibility y", "Connected Susceptibility y", "Binder Ratio y",
                        "|Magnetization|", "|Magnetization|^2", "|Magnetization|^4",
                        "Susceptibility", "Connected Susceptibility", "Binder Ratio",
                       ]
function loaddata(filename, obsnames)
    Ts = zeros(0)
    res = Parameter(n=>zeros(0) for n in obsnames)
    for line in eachline(filename)
        words = split(line)
        push!(Ts, parse(Float64, words[1]))
        for (i,n) in enumerate(obsnames)
            push!(res[n], parse(Float64, words[i+1]))
        end
    end
    return Ts, res
end

function parse_filename(filename, ::Union{Type{Ising}, Type{XY}})
    m = match(r"^J_([\d.-]*)__N_([\d.-]*).dat$", filename)
    if m == nothing
        return nothing
    end
    p = Parameter()
    p["J"] = parse(Float64, m.captures[1])
    p["N"] = parse(Int, m.captures[2])
    return p
end

function parse_filename(filename, ::Union{Type{Potts}, Type{Clock}})
    m = match(r"^Q_(\d*)__J_([\d.-]*)__N_(\d*).dat$", filename)
    if m == nothing
        return nothing
    end
    p = Parameter()
    p["Q"] = parse(Int, m.captures[1])
    p["J"] = parse(Float64, m.captures[2])
    p["N"] = parse(Int, m.captures[3])
    return p
end

@testset "$modelstr" for (modelstr, pnames, obsnames) in [("Ising", ("J", "N"), obsnames_ising),
                                                          ("Potts", ("Q", "J", "N"), obsnames_ising),
                                                          ("Clock", ("Q", "J", "N"), obsnames_clock),
                                                          ("XY", ("J", "N"), obsnames_clock),
                                                         ]
    model = eval(Symbol(modelstr))
    for filename in readdir(joinpath("ref", modelstr))
        p = parse_filename(filename, model)
        if p == nothing
            continue
        end
        testname = ""
        for pname in pnames
            testname *= "$(pname)=$(p[pname]) "
        end
        @testset "$testname" begin
            Ts, exacts = loaddata(joinpath("ref", modelstr, filename),obsnames)
            nT = length(Ts)
            p["Model"] = model
            p["Lattice"] = "fully connected graph"
            p["J"] = p["J"] / p["N"]
            p["MCS"] = MCS
            p["Thermalization"] = Therm
            @testset "$(upstr)" for upstr in ("local_update!",
                                              "SW_update!",
                                              "Wolff_update!",
                                             )
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
                        if !(p_value(mc1, exact) > alpha/nT || p_value(mc2, exact) > alpha/nT)
                            @show T
                            @show exact
                            @show mc1, p_value(mc1, exact)
                            @show mc2, p_value(mc2, exact)
                        end
                        @test p_value(mc1, exact) > alpha/nT || p_value(mc2, exact) > alpha/nT
                    end
                end
            end
        end
    end
end
