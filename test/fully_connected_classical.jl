function loaddata(filename)
    Ts = zeros(0)
    Es = zeros(0)
    for line in eachline(filename)
        words = split(line)
        push!(Ts, parse(words[1]))
        push!(Es, parse(words[2]))
    end
    return Ts, Es
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

@testset "$modelstr" for (modelstr, pnames) in [("Ising", ("J", "N")),
                                                ("Potts", ("Q", "J", "N")),
                                                ("Clock", ("Q", "J", "N")),
                                                ("XY", ("J", "N")),
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
            Ts, exacts = loaddata(joinpath("ref", modelstr, filename))
            nT = length(Ts)
            p["Model"] = model
            p["Lattice"] = fully_connected_lattice
            p["J"] = p["J"] / p["N"]
            p["MCS"] = MCS
            p["Thermalization"] = Therm
            @testset "UpdateMethod=$(up)" for up in (local_update!,
                                                     SW_update!,
                                                     Wolff_update!
                                                    )
                p["UpdateMethod"] = up
                for (T,exact) in zip(Ts,exacts)
                    srand(SEED)
                    p["T"] = T
                    res = runMC(p)
                    ene = res["Energy"]
                    if !(p_value(ene, exact) > alpha/nT)
                        ## Perform one more test since single MC test may fail.
                        res = runMC(p)
                        ene = res["Energy"]
                    end
                    if !(p_value(ene, exact) > alpha/nT)
                        @show T
                        @show exact
                        @show ene
                    end
                    @test p_value(ene, exact) > alpha/nT
                end
            end
        end
    end
end
