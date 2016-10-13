using SpinMonteCarlo

function gen_snapshot!(model::Model, T::Real; MCS::Integer=8192)
    for m in 1:MCS
        SW_update!(model, T)
    end
    return copy(model.spins)
end

function gen_snapshot!(model::Model, T::Real, N::Integer; MCS::Integer=8192)
    nsites = numsites(model.lat)
    X = zeros(nsites, N)
    for n in 1:N
        spins[:,n] = snapshot(model, T, MCS=MCS)
    end
    return X
end

function gensave_snapshot!(io::IO, model::Model, T::real, N::Integer=1; MCS::Integer=8192, sep::String=" ")
    nsites = numsites(model.lat)
    for n in 1:N
        for m in 1:MCS
            SW_update!(model, T)
        end
        print(io, spins[1])
        for s in 2:nsites
            print(io, sep, spins[s])
        end
        println(io)
    end
end

function gensave_snapshot!(filename::String, model::Model, T::real, N::Integer=1; MCS::Integer=8192, sep::String=" ", append::Bool=false)
    c = ifelse(append, "a", "w")
    open(filename, c) do io
        gensave_snapshot!(io, model, T, MCS=MCS, sep=sep)
    end
end

function load_snapshot_impl(io::IO, splitter)
    for line in eachline(io)
        if ismatch(r"^\s*($|#)", line)
            continue
        end
        produce(float(split(line, sep)))
    end
end
function load_snapshot_impl(filename::String, splitter)
    io = open(filename)
    for line in eachline(io)
        if ismatch(r"^\s*($|#)", line)
            continue
        end
        produce(float(split(line, sep)))
    end
    close(io)
end
const _default_delims = [' ','\t','\n','\v','\f','\r']
load_snapshot(io::IO, splitter = _default_delims) = @task load_snapshot_impl(io, splitter)
load_snapshot(filename::String, splitter = _default_delims) = @task load_snapshot_impl(filename, splitter)

