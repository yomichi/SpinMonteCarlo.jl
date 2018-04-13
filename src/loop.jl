"""
    loop_update!(model, param; measure::Bool=true)
    loop_update!(model, T::Real, Jz::Union{Real, AbstractArray}, Jxy::Union{Real, AbstractArray}, Gamma::Union{Real, AbstractArray}; measure::Bool=true)

updates spin configuration by loop algorithm 
under the temperature `T = param["T"]` and coupling constants `Jz, Jxy` and transverse field `Gamma`
"""
@inline function loop_update!(model::QuantumXXZ, param; measure::Bool=true)
    T = param["T"]
    if "J" in keys(param)
        Jz = Jxy = param["J"]
    else
        Jz = param["Jz"]
        Jxy = param["Jxy"]
    end
    Gamma = get(param, "Gamma", 0.0)
    return loop_update!(model, T, Jz, Jxy, Gamma, measure=measure)
end

@inline function loop_update!(model::QuantumXXZ, T::Real, Jz, Jxy, Gamma; measure::Bool=true)
    Jzs  = isa(Jz, Real) ? Jz  .* ones(numbondtypes(model)) : Jz
    Jxys  = isa(Jxy, Real) ? Jxy  .* ones(numbondtypes(model)) : Jxy
    Gammas  = isa(Gamma, Real) ? Gamma  .* ones(numsitetypes(model)) : Gamma
    return loop_update!(model, T, Jzs, Jxys, Gammas, measure=measure)
end
@inline function loop_update!(model::QuantumXXZ, T::Real, Jz, Jxy; measure::Bool=true)
    return loop_update!(model, T, Jz, Jxy, 0.0, measure=measure)
end
@inline function loop_update!(model::QuantumXXZ, T::Real, J; measure::Bool=true)
    return loop_update!(model, T, J, J, 0.0, measure=measure)
end

function loop_update!(model::QuantumXXZ, T::Real,
                      Jzs::AbstractArray, Jxys::AbstractArray, Gs::AbstractArray
                      ;
                      measure::Bool=true)
    rng = model.rng
    lo_types = [LO_FMLink, LO_AFLink, LO_Vertex, LO_Cross]
    nsites = numsites(model)
    S2 = model.S2
    nspins = nsites*S2
    nbonds = numbonds(model)
    nst = numsitetypes(model)
    nbt = numbondtypes(model)
    weights = zeros(4*nbt+nst)
    for i in 1:nbt
        nb = numbonds(model,i)*S2*S2
        z = nb*Jzs[i]
        x = nb*abs(Jxys[i])
        if z > x
            ## AntiFerroIsing like
            weights[(4i-3):(4i)] .= 0.5*[0.0, z-x, x, 0.0]
        elseif z < -x
            ## FerroIsing like
            weights[(4i-3):(4i)] .= 0.5*[-z-x, 0.0, 0.0, x]
        else
            ## XY like
            weights[(4i-3):(4i)] .= 0.25*[0.0, 0.0, x+z, x-z]
        end
    end
    for i in 1:nst
        ns = numsites(model,i)*S2
        weights[4nbt+i] = 0.5*ns*Gs[i]
    end
    accumulated_weights = cumsum(weights)
    op_dt = T/accumulated_weights[end]

    spins = model.spins[:]
    currents = collect(1:nspins)
    uf = UnionFind(nspins)

    ops = LocalOperator[]
    iops = 1 # index of the next operator in the original string
    t = randexp(rng)*op_dt
    while t <= 1.0 || iops <= length(model.ops)
        if iops > length(model.ops) || t < model.ops[iops].time
            ## INSERT
            ot = searchsortedfirst(accumulated_weights, rand(rng)*accumulated_weights[end])
            if ot <= 4nbt
                ## Bond
                bt = ceil(Int, ot/4)
                ot = mod1(ot,4)
                lo_type = lo_types[ot]
                b = bonds(model,bt)[rand(rng, 1:numbonds(model,bt))]
                s1 = source(model, b)
                s2 = target(model, b)
                ss1, ss2 = rand(rng, 1:S2, 2)
                if ifelse(spins[site2subspin(s1,ss1,S2)] == spins[site2subspin(s2,ss2,S2)],
                          lo_type == LO_FMLink || lo_type == LO_Cross,
                          lo_type == LO_AFLink || lo_type == LO_Vertex)
                    push!(ops, LocalOperator(lo_type, t, bond2subbond(b,ss1,ss2,S2)))
                    t += randexp(rng)*op_dt
                else
                    t += randexp(rng)*op_dt
                    continue
                end
            else
                ## Site
                st = ot-4nbt
                s = sites(model,st)[rand(rng, 1:numsites(model,st))]
                ss = rand(rng, 1:S2)
                push!(ops, LocalOperator(LO_Cut, t, site2subspin(s,ss,S2)))
                t += randexp(rng)*op_dt
            end
        else 
            op = model.ops[iops]
            iops += 1
            if op.isdiagonal
                ## REMOVE
                continue
            elseif op.op_type == LO_Cut
                push!(ops, op)
            else
                push!(ops, op)
                op = ops[end]
                subbond = op.space
                b,ss1,ss2 = subbond2bond(subbond,S2)
                bt = bondtype(model,b)
                ot = ifelse(rand(rng)*(weights[4bt-1]+weights[4bt])<weights[4bt-1], LO_Vertex, LO_Cross)
                op.op_type = ot
            end
        end
        
        op = ops[end]
        if op.op_type == LO_Cut
            subspin = op.space
            s,ss = subspin2site(subspin,S2)
            op.bottom_id = currents[subspin]
            c = addnode!(uf)
            currents[subspin] = c
            op.top_id = c
            spins[subspin] *= ifelse(op.isdiagonal, 1, -1)
        else
            subbond = op.space
            b,ss1,ss2 = subbond2bond(subbond,S2)
            s1 = source(model, b)
            s2 = target(model, b)
            subspin1 = site2subspin(s1,ss1,S2)
            subspin2 = site2subspin(s2,ss2,S2)
            if op.op_type == LO_FMLink || op.op_type == LO_AFLink
                c = unify!(uf, currents[subspin1], currents[subspin2])
                currents[subspin1] = currents[subspin2] = c
                op.bottom_id = op.top_id = c
            elseif op.op_type == LO_Cross
                op.bottom_id = currents[subspin1]
                op.top_id = currents[subspin2]
                spins[subspin1], spins[subspin2] = spins[subspin2], spins[subspin1]
                currents[subspin1], currents[subspin2] = currents[subspin2], currents[subspin1]
            else # if op.op_type == LO_Vertex
                unify!(uf, currents[subspin1], currents[subspin2])
                op.bottom_id = currents[subspin1]
                c = addnode!(uf)
                op.top_id = c
                currents[subspin1] = currents[subspin2] = c
                spins[subspin1] *= ifelse(op.isdiagonal, 1, -1)
                spins[subspin2] *= ifelse(op.isdiagonal, 1, -1)
            end
        end
    end # of while loop

    ## PBC for imaginary time axis
    subspin = 0
    for s in 1:nsites
        ups0 = zeros(Int,0)
        ups1 = zeros(Int,0)
        downs0 = zeros(Int,0)
        downs1 = zeros(Int,0)
        for ss in 1:S2
            subspin += 1
            push!(ifelse(spins[subspin]==1, ups0, downs0), subspin)
            push!(ifelse(model.spins[subspin]==1, ups1, downs1), subspin)
        end
        @assert length(ups0) == length(ups1)
        @assert length(downs0) == length(downs1)
        for (u, u2) in zip(ups0, shuffle(rng,ups1))
            unify!(uf, u, currents[u2])
        end
        for (d, d2) in zip(downs0, shuffle(rng,downs1))
            unify!(uf, d, currents[d2])
        end
    end

    model.ops = ops

    nc = clusterize!(uf)
    flips = rand(rng, [1,-1], nc)
    for s in 1:nspins
        model.spins[s] *= flips[clusterid(uf, s)]
    end
    for op in model.ops
        if op.op_type == LO_Cross || op.op_type == LO_Vertex || op.op_type == LO_Cut
            bid = clusterid(uf, op.bottom_id)
            tid = clusterid(uf, op.top_id)
            op.isdiagonal ⊻= (flips[bid] != flips[tid])
        end
    end

    res = Measurement()
    if measure
        M, M2, M4, E, E2, sgn = improved_estimate(model, T, Jzs, Jxys, Gs, uf)
        res["M"] = M
        res["M2"] = M2
        res["M4"] = M4
        res["E"] = E
        res["E2"] = E2
        res["Sign"] = sgn
    end
    return res
end
