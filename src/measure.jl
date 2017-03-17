"""
    measure(model::Ising, T::Real)
    measure(model::Potts, T::Real)

return magnetization density `M` and total energy `E`.
"""
function measure(model::Ising, T::Real)
    lat = model.lat
    nsites = numsites(lat)
    nbonds = numbonds(lat)

    M = mean(model.spins)
    E = 0.0
    @inbounds for b in 1:nbonds
        s1, s2 = source(lat, b), target(lat, b)
        E += ifelse(model.spins[s1] == model.spins[s2], -1.0, 1.0)
    end
    return M,E
end

function measure(model::Potts, T::Real)
    lat = model.lat
    nsites = numsites(lat)
    nbonds = numbonds(lat)
    invQ = 1.0/model.Q

    M = 0.0
    @inbounds for s in 1:nsites
        M += ifelse(model.spins[s]==1, 1.0-invQ, -invQ)
    end
    M /= nsites
    E = 0.0
    @inbounds for b in 1:nbonds
        s1, s2 = source(lat, b), target(lat, b)
        E -= ifelse(model.spins[s1] == model.spins[s2], 1.0, 0.0)
    end
    return M,E
end


"""
    measure(model::Clock, T::Real)
    measure(model::XY, T::Real)

return magnetization density `M`, total energy `E`, and helicity modulus `U`.
"""
function measure(model::Clock, T::Real)
    lat = model.lat
    nsites = numsites(lat)
    nbonds = numbonds(lat)
    D = dim(lat)
    invN = 1.0/nsites
    beta = 1.0/T

    M = zeros(2)
    E = 0.0
    U1 = zeros(2)
    U2 = zeros(2)

    @inbounds for s in 1:nsites
        M[1] += model.cosines[model.spins[s]]
        M[2] += model.sines[model.spins[s]]
    end

    @inbounds for b in 1:nbonds
        i, j = source(lat, b), target(lat,b)
        dir = bonddirectory(lat, b)
        dt = mod1(model.spins[j] - model.spins[i], model.Q)
        E -= model.cosines[dt]

        for d in 1:D
            U1[d] += model.cosines[dt] * dir[d]^2
            U2[d] += model.sines[dt] * dir[d]
        end
    end
    for d in 1:D
        M[d] *= invN
        U1[d] -= beta * U2[d]^2
        U1[d] *= invN
    end
    return M, E, U1
end

function measure(model::XY, T::Real)
    lat = model.lat
    nsites = numsites(lat)
    nbonds = numbonds(lat)
    D = dim(lat)
    invN = 1.0/nsites
    beta = 1.0/T

    M = zeros(2)
    E = 0.0
    U1 = zeros(2)
    U2 = zeros(2)

    @inbounds for s in 1:nsites
        M[1] += cospi(2model.spins[s])
        M[2] += sinpi(2model.spins[s])
    end

    @inbounds for b in 1:nbonds
        i, j = source(lat, b), target(lat,b)
        dir = bonddirectory(lat, b)
        dt = mod(model.spins[j] - model.spins[i] + 2.0, 1.0)
        dt = ifelse(0.5 <= dt, dt - 1.0, dt)
        E -= cospi(2dt)

        for d in 1:D
            U1[d] += cospi(2dt) * dir[d]^2
            U2[d] += sinpi(2dt) * dir[d]
        end
    end
    for d in 1:D
        M[d] *= invN
        U1[d] -= beta * U2[d]^2
        U1[d] *= invN
    end
    return M, E, U1
end

function measure(model::QuantumLocalZ2Model, uf::UnionFind)
    nsites = numsites(model.lat)
    nbonds = numbonds(model.lat)
    nc = numclusters(uf)

    spins = model.spins[:]
    ms = zeros(nc)
    @inbounds for op in model.ops
        if op.op_type == LO_Cut
            s = op.space
            bid = clusterid(uf, op.bottom_id)
            tid = clusterid(uf, op.top_id)
            ms[bid] += spins[s]*op.time
            ms[tid] -= spins[s]*op.time
            spins[s] *= ifelse(op.isdiagonal, 1, -1)
        elseif op.op_type == LO_Vertex
            b = op.space
            s1 = source(model.lat, b)
            s2 = target(model.lat, b)
            bid = clusterid(uf, op.bottom_id)
            tid = clusterid(uf, op.top_id)
            ms[bid] += (spins[s1]+spins[s2])*op.time
            ms[tid] -= (spins[s1]+spins[s2])*op.time
            spins[s1] *= ifelse(op.isdiagonal, 1, -1)
            spins[s2] *= ifelse(op.isdiagonal, 1, -1)
        elseif op.op_type == LO_Cross
            b = op.space
            s1 = source(model.lat, b)
            s2 = target(model.lat, b)
            bid = clusterid(uf, op.bottom_id)
            tid = clusterid(uf, op.top_id)
            ms[bid] += (spins[s1]-spins[s2])*op.time
            ms[tid] -= (spins[s1]-spins[s2])*op.time
            spins[s1] *= ifelse(op.isdiagonal, 1, -1)
            spins[s2] *= ifelse(op.isdiagonal, 1, -1)
        end
    end
    for s in 1:nsites
        i = clusterid(uf, s)
        ms[i] += model.spins[s]
    end

    coeff = 0.5/nsites
    @inbounds for i in 1:nc
        ms[i] *= coeff
    end

    M = 0.0
    M2 = 0.0
    M4 = 0.0
    for m in ms
        M += m
        m2 = m*m
        M4 += m2*m2 + 6M2*m2
        M2 += m2
    end
    return M, M2, M4
end
