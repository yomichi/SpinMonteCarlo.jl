@doc """
    improved_estimator(model::QuantumXXZ, T::Real, Js::AbstractArray, uf::UnionFind)

Returns the following observables as `Dict{String, Any}` using loop information `uf`

# Observables
- `"Sign"`
    - Sign of the weight function
- `"Sign * Energy"`
    - Energy per spin (site)
- `"Sign * Energy^2"`
- `"Sign * Magnetization"`
    - Total magnetization (Sz) per spin (site)
- `"Sign * |Magnetization|"`
- `"Sign * Magnetization^2"`
- `"Sign * Magnetization^4"`
"""
function improved_estimator(model::QuantumXXZ, T::Real, Jzs::AbstractArray,
                            Jxys::AbstractArray, Gs::AbstractArray, uf::UnionFind)
    S2 = model.S2
    nsites = numsites(model)
    nspins = nsites * S2
    nbonds = numbonds(model)
    nc = numclusters(uf)

    ms = zeros(nc)
    for s in 1:nspins
        i = clusterid(uf, s)
        ms[i] += model.spins[s]
    end
    ms .*= 0.5 / nsites

    E0 = 0.0
    for st in 1:numsitetypes(model)
        E0 += 0.5 * numsites(model, st) * Gs[st] * S2
    end
    for bt in 1:numbondtypes(model)
        nb = numbonds(model, bt) * S2 * S2
        z = nb * Jzs[bt]
        x = nb * abs(Jxys[bt])
        if z > x
            ## AntiFerroIsing like
            E0 += 0.25z
        elseif z < -x
            ## FerroIsing like
            E0 -= 0.25z
        else
            ## XY like
            E0 += 0.25x
        end
    end
    nops = length(model.ops)
    E = E0 - nops * T
    E2 = nops * (nops - 1) * T^2 - 2 * E0 * T * nops + E0^2
    E /= nsites
    E2 /= nsites^2

    M = 0.0
    M2 = 0.0
    M4 = 0.0
    for m in ms
        M += m
        m2 = m * m
        M4 += m2 * m2 + 6M2 * m2
        M2 += m2
    end

    sgn = 1.0
    for op in model.ops
        if !op.isdiagonal
            if op.let_type == LET_Vertex || op.let_type == LET_Cross
                b = op.space - nsites
                bt = bondtype(model, b)
                sgn *= ifelse(Jxys[bt] > 0.0, -1.0, 1.0)
            end
        end
    end

    res = Measurement()
    res["Sign * Magnetization"] = sgn * M
    res["Sign * |Magnetization|"] = sgn * abs(M)
    res["Sign * Magnetization^2"] = sgn * M2
    res["Sign * Magnetization^4"] = sgn * M4
    res["Sign * Energy"] = sgn * E
    res["Sign * Energy^2"] = sgn * E2
    res["Sign"] = sgn

    return res
end

default_estimator(model::QuantumXXZ, update) = improved_estimator
