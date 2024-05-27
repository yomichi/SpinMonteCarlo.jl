export MCObservableSet

"""
    MCObservableSet{Obs}

    Alias of `Dict{String, Obs}` where `Obs` is a subtype of `MCObservable`.
"""
const MCObservableSet{Obs<:MCObservable} = Dict{String,Obs}

export makeMCObservable!

"""
    makeMCObservable!(oset::MCObservableSet{Obs}, name::String)

    Create an observable with the name `name` in the observable set `oset`.
"""
function makeMCObservable!(oset::MCObservableSet{Obs},
                           name::String) where {Obs<:MCObservable}
    if haskey(oset, name)
        warn("""Observable "$name" already exists. (Skipped)""")
    else
        oset[name] = Obs()
    end
end

function reset!(oset::MCObservableSet)
    for v in values(oset)
        reset!(oset)
    end
end

function show(io::IO, obs::MCObservableSet, sorted::Bool=true)
    if sorted
        ks = sort([k for k in keys(obs)])
        for k in ks
            println(io, k, " : ", obs[k])
        end
    else
        for k in keys(obs)
            println(io, k, " : ", obs[k])
        end
    end
end
