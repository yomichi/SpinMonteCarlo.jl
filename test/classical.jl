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
                        "SUsceptibility", "Connected Susceptibility", "Binder Ratio",
                       ]
function loaddata(filename, obsnames)
    Ts = zeros(0)
    res = Dict(n=>zeros(0) for n in obsnames)
    for line in eachline(filename)
        words = split(line)
        push!(Ts, parse(words[1]))
        for (i,n) in enumerate(obsnames)
            push!(res[n], parse(words[i+1]))
        end
    end
    return Ts, res
end

function parse_filename(filename, ::Union{Type{Ising}, Type{XY}})
    m = match(r"^J_([\d.-]*)__N_([\d.-]*).dat$", filename)
    if m == nothing
        return nothing
    end
    p = Dict()
    p["J"] = parse(m.captures[1])
    p["N"] = parse(Int, m.captures[2])
    return p
end

function parse_filename(filename, ::Union{Type{Potts}, Type{Clock}})
    m = match(r"^Q_(\d*)__J_([\d.-]*)__N_(\d*).dat$", filename)
    if m == nothing
        return nothing
    end
    p = Dict()
    p["Q"] = parse(Int, m.captures[1])
    p["J"] = parse(m.captures[2])
    p["N"] = parse(Int, m.captures[3])
    return p
end

@testset "$modelstr" for (modelstr, pnames, obsnames) in [("Ising", ("J", "N"), obsnames_ising),
                                                          ("Potts", ("Q", "J", "N"), obsnames_ising),
                                                          # ("Clock", ("Q", "J", "N"), obsnames_clock),
                                                          # ("XY", ("J", "N"), obsnames_clock),
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
            p["Lattice"] = fully_connected_lattice
            p["J"] = p["J"] / p["N"]
            p["MCS"] = MCS
            p["Thermalization"] = Therm
            @testset "UpdateMethod=$(upstr)" for upstr in ("local_update!",
                                                           "SW_update!",
                                                           "Wolff_update!",
                                                          )
                p["UpdateMethod"] = eval(Symbol(upstr))
                @testset "$n" for n in obsnames
                    for (T,exact) in zip(Ts,exacts[n])
                        srand(SEED)
                        p["T"] = T
                        res1 = runMC(p)
                        res2 = runMC(p)
                        ## single MC test may fail.
                        mc1 = res1[n]
                        mc2 = res2[n]
                        ex = exact
                        if !(p_value(mc1, exact) > alpha/nT || p_value(mc2, exact) > alpha/nT)
                            @show T
                            @show exact
                            @show mc1
                            @show mc2
                        end
                        @test p_value(mc1, exact) > alpha/nT || p_value(mc2, exact) > alpha/nT
                    end
                end
            end
        end
    end
end
