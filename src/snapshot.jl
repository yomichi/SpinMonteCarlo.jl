@doc """
    gen_snapshot!(model, T, [N=1])

generate and return `N` snapshots (spin configuration).

# Arguments
* `model::Model` : model
* `T::Real` : temperature
* `N::Integer=1` : the number of snapshots to be generated
* `MCS::Integer=8192` : the number of Monte Carlo steps followed by snapshot generation
"""
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
        spins[:,n] = gen_snapshot!(model, T, MCS=MCS)
    end
    return X
end

@doc """
    gensave_snapshot!(io, model, T, [N=1])
    gensave_snapshot!(filename, model, T, [N=1])

generate and write `N` snapshots into `io` or `filename`.

# Arguments
* `io::IO` : output stream where snapshots will be written
* `filename::String` : the name of file where snapshots will be written
* `model::Model` : model to be simulated
* `T::Real` : temperature
* `N::Integer=1` : the number of snapshots to be generated
* `MCS::Integer=1` : the number of Monte Carlo steps followed by snapshot generation
* `sep::String=" "` : field separator of snapshot
"""
function gensave_snapshot!(io::IO, model::Model, T::Real, N::Integer=1; MCS::Integer=8192, sep::String=" ")
    nsites = numsites(model.lat)
    for n in 1:N
        for m in 1:MCS
            SW_update!(model, T)
        end
        spins = model.spins
        print(io, spins[1])
        for s in 2:nsites
            print(io, sep, spins[s])
        end
        println(io)
    end
end

function gensave_snapshot!(filename::String, model::Model, T::Real, N::Integer=1; MCS::Integer=8192, sep::String=" ", append::Bool=false)
    c = ifelse(append, "a", "w")
    open(filename, c) do io
        gensave_snapshot!(io, model, T, N, MCS=MCS, sep=sep)
    end
end

function load_snapshot_impl(io::IO, splitter)
    for line in eachline(io)
        if ismatch(r"^\s*($|#)", line)
            continue
        end
        produce(float(split(strip(line), splitter)))
    end
end
function load_snapshot_impl(filename::String, splitter)
    io = open(filename)
    for line in eachline(io)
        if ismatch(r"^\s*($|#)", line)
            continue
        end
        produce(float(split(strip(line), splitter)))
    end
    close(io)
end
const _default_delims = [' ','\t','\n','\v','\f','\r']

@doc """
    load_snapshot(io, splitter)
    load_snapshot(filename, splitter)

load snapshot from `io` or `filename`.
`splitter` is the set of field splitter and will be passed `Base.split` as the second argument (see the doc of `Base.split`).
"""
load_snapshot(io::IO, splitter = _default_delims) = @task load_snapshot_impl(io, splitter)
load_snapshot(filename::String, splitter = _default_delims) = @task load_snapshot_impl(filename, splitter)

