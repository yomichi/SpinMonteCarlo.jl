function loop_update!(model::TransverseFieldIsing, T::Real, J::Real, gamma::Real; measure::Bool=true)
    nsites = numsites(model.lat)
    nbonds = numbonds(model.lat)
    ising_weight = 0.5*abs(J)*nbonds
    field_weight = 0.5*abs(gamma)*nsites
    op_dt = T/(ising_weight + field_weight)
    r_ising = ising_weight/(ising_weight+field_weight)

    currents = collect(1:nsites)
    uf = UnionFind(nsites)

    ops = LocalOperator[]
    iops = 1 # index of the next operator
    t = randexp()*op_dt
    while t <= 1.0 || iops <= length(model.ops)
        if iops > length(model.ops) || t < model.ops[iops].time
            # INSERT

            if rand() < r_ising
                b = rand(1:nbonds)
                s1 = source(model.lat, b)
                s2 = target(model.lat, b)
                if J * model.spins[s1] * model.spins[s2] < 0.0
                    push!(ops, LocalOperator(LO_Link, t, b))
                    t += randexp()*op_dt
                else
                    t += randexp()*op_dt
                    continue
                end
            else
                s = rand(1:nsites)
                push!(ops, LocalOperator(LO_Cut, t, s))
                t += randexp()*op_dt
            end
        else 
            if ! model.ops[iops].isdiagonal
                push!(ops, model.ops[iops])
                iops += 1
            else
                iops += 1
                continue
            end
        end
        
        op = ops[end]
        if op.op_type == LO_Link
            b = op.space
            s1 = source(model.lat, b)
            s2 = target(model.lat, b)
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
    end
    for s in 1:nsites
        unify!(uf, s, currents[s])
    end

    model.ops = ops

    nc = clusterize!(uf)
    flips = rand([1,-1], nc)
    for s in 1:nsites
        model.spins[s] = flips[clusterid(uf, s)]
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
        M, M2, M4 = improved_estimate(model, uf)
        res[:M] = M
        res[:M2] = M2
        res[:M4] = M4
    end
    return res
end
