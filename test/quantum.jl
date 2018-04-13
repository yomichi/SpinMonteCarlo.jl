const obsnames_XXZ = ["Energy", "Energy^2", "Specific Heat",
                      "Magnetization^2", "Magnetization^4",
                      "Binder Ratio", "Susceptibility",
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

function parse_filename(filename)
    m = match(r"^S_([\d.-]*)__Jz_([\d.-]*)__Jxy_([\d.-]*)__G_([\d.-]*)__L_([\d.-]*).dat$", filename)
    if m == nothing
        return nothing
    end
    p = Dict()
    p["S"] = parse(m.captures[1])
    p["Jz"] = parse(m.captures[2])
    p["Jxy"] = parse(m.captures[3])
    p["Gamma"] = parse(m.captures[4])
    p["L"] = parse(Int, m.captures[5])
    return p
end

@testset "$modelstr" for (modelstr, pnames, obsnames, updatestrs) in [("QuantumXXZ",
                                                                       ("S", "Jz", "Jxy", "Gamma", "L"),
                                                                       obsnames_XXZ,
                                                                       ("loop_update!",),
                                                                      ),
                                                                     ]
    model = eval(Symbol(modelstr))
    @testset "$(latticestr)" for latticestr in ["chain_lattice"]
        lattice = eval(Symbol(latticestr))
        for filename in readdir(joinpath("ref", modelstr, latticestr))
            p = parse_filename(filename)
            if p == nothing
                continue
            end
            testname = ""
            for pname in pnames
                testname *= "$(pname)=$(p[pname]) "
            end
            @testset "$testname" begin
                Ts, exacts = loaddata(joinpath("ref", modelstr, latticestr, filename),obsnames)
                nT = length(Ts)
                p["Model"] = model
                p["Lattice"] = lattice
                p["MCS"] = MCS
                p["Thermalization"] = Therm
                @testset "$updatestr" for updatestr in updatestrs
                    p["Update Method"] = eval(Symbol(updatestr))
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
                            if  mean(r1["Sign"]) < 1.0 && (!(isfinite(mean(mc1))) || !(isfinite(mean(mc2))))
                                ## Sign problem makes test very difficult...
                                continue
                            end
                            # @show exact, mc1, mc2
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
end
