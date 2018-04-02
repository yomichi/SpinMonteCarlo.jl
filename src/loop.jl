function loop_update!(model::TransverseFieldIsing, T::Real, J::Real, G::Real; measure::Bool=true)
    Js = J * ones(numbondtypes(model))
    Gs = G * ones(numsitetypes(model))
    return loop_update!(model, T, Js, Gs, measure=measure)
end
function loop_update!(model::TransverseFieldIsing, T::Real, J::Real, Gs::AbstractArray; measure::Bool=true)
    Js = J * ones(numbondtypes(model))
    return loop_update!(model, T, Js, Gs, measure=measure)
end
function loop_update!(model::TransverseFieldIsing, T::Real, Js::AbstractArray, G::Real; measure::Bool=true)
    Gs = G * ones(numsitetypes(model))
    return loop_update!(model, T, Js, Gs, measure=measure)
end
function loop_update!(model::TransverseFieldIsing, T::Real, Js::AbstractArray, Gs::AbstractArray; measure::Bool=true)
    nsites = numsites(model)
    nbonds = numbonds(model)

    nbt = numbondtypes(model)
    ising_weights = zeros(nbt)
    for i in 1:nbt
        ising_weights[i] = 0.5*abs(Js[i])*numbonds(model,i)
    end
    cumsum!(ising_weights, ising_weights)

    nst = numsitetypes(model)
    field_weights = zeros(nst)
    for i in 1:nst
        field_weights[i] = 0.5*abs(Gs[i])*numsites(model,i)
    end
    cumsum!(field_weights, field_weights)

    op_dt = T/(ising_weights[end] + field_weights[end])
    r_ising = ising_weights[end]/(ising_weights[end]+field_weights[end])

    currents = collect(1:nsites)
    uf = UnionFind(nsites)

    ops = LocalOperator[]
    iops = 1 # index of the next operator
    t = randexp()*op_dt
    while t <= 1.0 || iops <= length(model.ops)
        if iops > length(model.ops) || t < model.ops[iops].time
            ## INSERT

            if rand() < r_ising
                r = rand()*ising_weights[end]
                bt = searchsortedfirst(ising_weights,r)
                b = bonds(model,bt)[rand(1:numbonds(model,bt))]
                s1 = source(model, b)
                s2 = target(model, b)
                if Js[bt] * model.spins[s1] * model.spins[s2] < 0.0
                    push!(ops, LocalOperator(ifelse(Js[bt] > 0.0, LO_AFLink, LO_FMLink), t, b))
                    t += randexp()*op_dt
                else
                    t += randexp()*op_dt
                    continue
                end
            else
                r = rand()*field_weights[end]
                st = searchsortedfirst(field_weights,r)
                s = sites(model,st)[rand(1:numsites(model,st))]
                push!(ops, LocalOperator(LO_Cut, t, s))
                t += randexp()*op_dt
            end
        else 
            ## REMOVE

            if ! model.ops[iops].isdiagonal
                push!(ops, model.ops[iops])
                iops += 1
            else
                iops += 1
                continue
            end
        end
        
        op = ops[end]
        if op.op_type == LO_FMLink || op.op_type == LO_AFLink
            b = op.space
            s1 = source(model, b)
            s2 = target(model, b)
            c = unify!(uf, currents[s1], currents[s2])
            currents[s1] = currents[s2] = c
            op.bottom_id = op.top_id = c
        elseif op.op_type == LO_Cut
            s = op.space
            op.bottom_id = currents[s]
            c = addnode!(uf)
            op.top_id = c
            currents[s] = c
            model.spins[s] *= ifelse(op.isdiagonal, 1, -1)
        end
    end # of while loop

    ## PBC for imaginary time axis
    for s in 1:nsites
        unify!(uf, s, currents[s])
    end

    model.ops = ops

    nc = clusterize!(uf)
    flips = rand([1,-1], nc)
    for s in 1:nsites
        model.spins[s] *= flips[clusterid(uf, s)]
    end
    for op in model.ops
        if op.op_type == LO_Cut
            bid = clusterid(uf, op.bottom_id)
            tid = clusterid(uf, op.top_id)
            op.isdiagonal = (flips[bid] == flips[tid])
        end
    end

    res = Measurement()
    if measure
        M, M2, M4, E, E2 = improved_estimate(model, T, Js, Gs, uf)
        res["M"] = M
        res["M2"] = M2
        res["M4"] = M4
        res["E"] = E
        res["E2"] = E2
    end
    return res
end

function loop_update!(model::QuantumXXZ, T::Real, Jz, Jxy, G; measure::Bool=true)
    Jzs  = isa(Jz, Real) ? Jz  .* ones(numbondtypes(model)) : Jz
    Jxys  = isa(Jxy, Real) ? Jxy  .* ones(numbondtypes(model)) : Jxy
    Gs  = isa(G, Real) ? G  .* ones(numsitetypes(model)) : G
    return loop_update!(model, T, Jzs, Jxys, Gs, measure=measure)
end
function loop_update!(model::QuantumXXZ, T::Real, Jz, Jxy; measure::Bool=true)
    return loop_update!(model, T, Jz, Jxy, 0.0, measure=measure)
end
function loop_update!(model::QuantumXXZ, T::Real, J; measure::Bool=true)
    return loop_update!(model, T, J, J, 0.0, measure=measure)
end

function loop_update!(model::QuantumXXZ, T::Real,
                      Jzs::AbstractArray, Jxys::AbstractArray, Gs::AbstractArray
                      ;
                      measure::Bool=true)
    lo_types = [LO_FMLink, LO_AFLink, LO_Vertex, LO_Cross]
    nsites = numsites(model)
    S2 = model.S2
    nspins = nsites*S2
    nbonds = numbonds(model)
    nst = numsitetypes(model)
    nbt = numbondtypes(model)
    shift = 0.0
    weights = zeros(4*nbt+nst)
    for i in 1:nbt
        nb = numbonds(model,i)*S2*S2
        z = nb*Jzs[i]
        x = nb*abs(Jxys[i])
        if z > x
            ## AntiFerroIsing like
            weights[(4i-3):(4i)] .= 0.5*[0.0, z-x, x, 0.0]
            shift += 0.25z
        elseif z < -x
            ## FerroIsing like
            weights[(4i-3):(4i)] .= 0.5*[-z-x, 0.0, 0.0, x]
            shift -= 0.25z
        else
            ## XY like
            weights[(4i-3):(4i)] .= 0.25*[0.0, 0.0, x+z, x-z]
            shift += 0.25x
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
    t = randexp()*op_dt
    while t <= 1.0 || iops <= length(model.ops)
        if iops > length(model.ops) || t < model.ops[iops].time
            ## INSERT
            ot = searchsortedfirst(accumulated_weights, rand()*accumulated_weights[end])
            if ot <= 4nbt
                ## Bond
                bt = ceil(Int, ot/4)
                ot = mod1(ot,4)
                lo_type = lo_types[ot]
                b = bonds(model,bt)[rand(1:numbonds(model,bt))]
                s1 = source(model, b)
                s2 = target(model, b)
                ss1, ss2 = rand(1:S2, 2)
                if ifelse(spins[site2subspin(s1,ss1,S2)] == spins[site2subspin(s2,ss2,S2)],
                          lo_type == LO_FMLink || lo_type == LO_Cross,
                          lo_type == LO_AFLink || lo_type == LO_Vertex)
                    push!(ops, LocalOperator(lo_type, t, bond2subbond(b,ss1,ss2,S2)))
                    t += randexp()*op_dt
                else
                    t += randexp()*op_dt
                    continue
                end
            else
                ## Site
                st = ot-4nbt
                s = sites(model,st)[rand(1:numsites(model,st))]
                ss = rand(1:S2)
                push!(ops, LocalOperator(LO_Cut, t, site2subspin(s,ss,S2)))
                t += randexp()*op_dt
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
                ot = ifelse(rand()*(weights[4bt-1]+weights[4bt])<weights[4bt-1], LO_Vertex, LO_Cross)
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
        for (u, u2) in zip(ups0, shuffle(ups1))
            unify!(uf, u, currents[u2])
        end
        for (d, d2) in zip(downs0, shuffle(downs1))
            unify!(uf, d, currents[d2])
        end
    end

    model.ops = ops

    nc = clusterize!(uf)
    flips = rand([1,-1], nc)
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
        M, M2, M4, E, E2 = improved_estimate(model, T, Jzs, Jxys, Gs, uf)
        res["M"] = M
        res["M2"] = M2
        res["M4"] = M4
        res["E"] = E
        res["E2"] = E2
    end
    return res
end
