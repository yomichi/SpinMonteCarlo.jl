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
