import Random: shuffle

@generated function myshuffle(xs)
    canShuffleEmpty = try
        shuffle([])
        true
    catch
        false
    end
    if canShuffleEmpty
        :(shuffle(xs))
    else
        :(length(xs)>0 ? shuffle(xs) : xs)
    end
end

@generated function myshuffle(rng, xs)
    canShuffleEmpty = try
        shuffle([])
        true
    catch
        false
    end
    if canShuffleEmpty
        :(shuffle(rng, xs))
    else
        :(length(xs)>0 ? shuffle(rng, xs) : xs)
    end
end

@doc """
    loop_update!(model, param::Parameter)
    loop_update!(model, T::Real,
                 Jz::AbstractArray,
                 Jxy::AbstractArray,
                 Gamma:AbstractArray)

Updates spin configuration by loop algorithm 
under the temperature `T = param["T"]` and coupling constants `Jz, Jxy` and transverse field `Gamma`
"""
@inline function loop_update!(model::QuantumXXZ, param::Parameter)
    p = convert_parameter(model, param)
    return loop_update!(model, p...)
end

function loop_update!(model::QuantumXXZ, T::Real,
                      Jzs::AbstractArray, Jxys::AbstractArray, Gs::AbstractArray)
    rng = model.rng
    lo_types = [LET_FMLink, LET_AFLink, LET_Vertex, LET_Cross]
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

    ops = LocalLoopOperator[]
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
                s1,s2 = source(b), target(b)
                ss1, ss2 = rand(rng, 1:S2, 2)
                if ifelse(spins[site2subspin(s1,ss1,S2)] == spins[site2subspin(s2,ss2,S2)],
                          lo_type == LET_FMLink || lo_type == LET_Cross,
                          lo_type == LET_AFLink || lo_type == LET_Vertex)
                    push!(ops, LocalLoopOperator(lo_type, t, nsites+b.id, (ss1, ss2)))
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
                push!(ops, LocalLoopOperator(LET_Cut, t, s, (ss,ss)))
                t += randexp(rng)*op_dt
            end
        else 
            op = model.ops[iops]
            iops += 1
            if op.isdiagonal
                ## REMOVE
                continue
            elseif op.let_type == LET_Cut
                push!(ops, op)
            else
                push!(ops, op)
                op = ops[end]
                b = op.space-nsites
                ss1,ss2 = op.subspace
                bt = bondtype(model,b)
                otype = ifelse(rand(rng)*(weights[4bt-1]+weights[4bt])<weights[4bt-1], LET_Vertex, LET_Cross)
                op.let_type = otype
            end
        end
        
        op = ops[end]
        if op.let_type == LET_Cut
            s = op.space
            ss = op.subspace[1]
            subspin = site2subspin(s,ss,S2)
            op.bottom_id = currents[subspin]
            c = addnode!(uf)
            currents[subspin] = c
            op.top_id = c
            spins[subspin] *= ifelse(op.isdiagonal, 1, -1)
        else
            b = op.space - nsites
            ss1,ss2 = op.subspace
            s1 = source(model, b)
            s2 = target(model, b)
            subspin1 = site2subspin(s1,ss1,S2)
            subspin2 = site2subspin(s2,ss2,S2)
            if op.let_type == LET_FMLink || op.let_type == LET_AFLink
                c = unify!(uf, currents[subspin1], currents[subspin2])
                currents[subspin1] = currents[subspin2] = c
                op.bottom_id = op.top_id = c
            elseif op.let_type == LET_Cross
                op.bottom_id = currents[subspin1]
                op.top_id = currents[subspin2]
                spins[subspin1], spins[subspin2] = spins[subspin2], spins[subspin1]
                currents[subspin1], currents[subspin2] = currents[subspin2], currents[subspin1]
            else # if op.let_type == LET_Vertex
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
        for (u, u2) in zip(ups0, myshuffle(rng,ups1))
            unify!(uf, u, currents[u2])
        end
        for (d, d2) in zip(downs0, myshuffle(rng,downs1))
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
        if op.let_type == LET_Cross || op.let_type == LET_Vertex || op.let_type == LET_Cut
            bid = clusterid(uf, op.bottom_id)
            tid = clusterid(uf, op.top_id)
            op.isdiagonal ⊻= (flips[bid] != flips[tid])
        end
    end

    return uf
end
