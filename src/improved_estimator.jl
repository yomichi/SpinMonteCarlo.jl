"""
    improved_estimate(model::Ising, T::Real, sw::SWInfo)

    return `M`, `M2`, `M4`, `E`, `E2`.

    `M`  : magnetization density `M`
    `M2` : square of mag density
    `M4` : biquadratic of mag density
    `E`  : energy density
    `E2` : square of energy density
"""
function improved_estimate(model::Ising, T::Real, sw::SWInfo)
    nsites = numsites(model)
    nbonds = numbonds(model)
    nc = numclusters(sw)
    invV = 1.0/nsites

    ## magnetization
    M = 0.0
    M2 = 0.0
    M4 = 0.0
    for (m,s) in zip(sw.clustersize, sw.clusterspin)
        M += m*invV*s
        m2 = (m*invV)^2
        M4 += m2*m2 + 6M2*m2
        M2 += m2
    end

    # energy
    beta = 1.0/T
    lambda = expm1(2beta)
    A = exp(2beta)/lambda
    n = sw.activated_bonds
    E = (-2A*n + nbonds)*invV
    E2 = (4A*A*n*(n-1) + 4A*n*(1-nbonds) + nbonds*nbonds) * invV*invV

    return M, M2, M4, E, E2
end

"""
    improved_estimate(model::Potts, T::Real, sw::SWInfo)

    return `M`, `M2`, `M4`, `E`, `E2`.

    `M`  : magnetization density `M`
    `M2` : square of mag density
    `M4` : biquadratic of mag density
    `E`  : energy density
    `E2` : square of energy density

    local magnetization `M_i` is defined by local spin variable `s_i` as `M_i = \\delta_{s_i, 1}-1/q`.
"""
function improved_estimate(model::Potts, T::Real, sw::SWInfo)
    nsites = numsites(model)
    Q = model.Q
    nc = numclusters(sw)
    invV = 1.0/nsites

    # magnetization
    I2 = (Q-1)/(Q*Q)
    I4 = (Q-1)*((Q-1)^3+1)/(Q^5)
    M = 0.0
    M2 = 0.0
    M4 = 0.0
    s2 = -1.0/Q
    s1 = 1.0+s2
    for (m,s) in zip(sw.clustersize, sw.clusterspin)
        M += (m*invV) * ifelse(s==1, s1, s2)
        m2 = (m*invV)^2
        M4 += I4*m2*m2 + 6*I2*M2*m2
        M2 += I2*m2
    end

    # energy
    beta = 1.0/T
    lambda = expm1(beta)
    A = exp(beta)/lambda
    n = sw.activated_bonds
    E = -A*n
    E2 = A*A*n*(n-1) + A*n
    E *= invV
    E2 *= invV*invV
    return M, M2, M4, E, E2
end

function improved_estimate(model::QuantumLocalZ2Model, uf::UnionFind)
    nsites = numsites(model)
    nbonds = numbonds(model)
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
            s1 = source(model, b)
            s2 = target(model, b)
            bid = clusterid(uf, op.bottom_id)
            tid = clusterid(uf, op.top_id)
            ms[bid] += (spins[s1]+spins[s2])*op.time
            ms[tid] -= (spins[s1]+spins[s2])*op.time
            spins[s1] *= ifelse(op.isdiagonal, 1, -1)
            spins[s2] *= ifelse(op.isdiagonal, 1, -1)
        elseif op.op_type == LO_Cross
            b = op.space
            s1 = source(model, b)
            s2 = target(model, b)
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
