function print_result(io::IO, params::Array, obs::Array)
    if isa(params[1]["Model"], Union{Type{Ising}, Type{XY}})
        names = ["L", "T"]
    elseif isa(params[1]["Model"], Union{Type{Potts}, Type{Clock}})
        names = ["Q", "L", "T"]
    elseif isa(params[1]["Model"], Union{Type{QuantumXXZ}})
        names = ["S2", "L", "T", "Jz", "Jxy", "Gamma"]
    end
    print_result(io, params, obs, names)
end

function print_result(io::IO, params, obs::Array, param_names::Array)
    i  = 1
    for name in param_names
        println(io, "# $i : $name")
        i += 1
    end
    ks = collect(keys(obs[1]))
    sort!(ks)
    for name in ks
        println(io, "# $i, $(i+1) : $name")
        i += 2
    end

    for (p, o) in zip(params, obs)
        for name in param_names
            print(io, p[name], " ")
        end

        for name in ks
            print(io, mean(o[name]), " ")
            print(io, stderror(o[name]), " ")
        end
        println(io)
    end
end

