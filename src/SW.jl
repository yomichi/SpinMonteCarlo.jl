type SWInfo
    activated_bonds :: Int
    clustersize :: Vector{Int}
    clusterspin :: Vector{Int}
end

numclusters(sw::SWInfo) = length(sw.clustersize)

"""
    magnetizations(sw::SWInfo, model)

return square and biquadratic of magnetization density.
"""
function magnetizations(sw::SWInfo, model::Ising)
    nsites = numsites(model.lat)
    nc = numclusters(sw)
    invV = 1.0/nsites
    M2 = 0.0
    M4 = 0.0
    for m in sw.clustersize
        m2 = (m*invV)^2
        M4 += m2*m2 + 6M2*m2
        M2 += m2
    end
    return M2, M4
end

function magnetizations(sw::SWInfo, model::Potts)
    nsites = numsites(model.lat)
    Q = model.Q
    nc = numclusters(sw)
    invV = 1.0/nsites
    I2 = (Q-1)/(Q*Q)
    I4 = (Q-1)*((Q-1)^3+1)/(Q^5)
    M2 = 0.0
    M4 = 0.0
    for m in sw.clustersize
        m2 = (m*invV)^2
        M4 += I4*m2*m2 + 6*M2*m2
        M2 += I2*m2
    end
    return M2, M4
end

"""
    SW_update!(model, T::Real)
    
update spin configuration by Swendsen-Wang algorithm under the temperature `T`, and return clusters.
"""
function SW_update!(model::Ising, T::Real)
    p = -expm1(-2.0/T)
    nsites = numsites(model.lat)
    nbonds = numbonds(model.lat)
    activated_bonds = 0
    uf = UnionFind(nsites)
    @inbounds for bond in 1:nbonds
        s1,s2 = source(model.lat, bond), target(model.lat, bond)
        if model.spins[s1] == model.spins[s2] && rand() < p
            activated_bonds += 1
            unify!(uf, s1,s2)
        end
    end
    nc = clusterize!(uf)
    clustersize = ones(Int, nc)
    clusterspin = rand([1,-1], nc)

    @inbounds for site in 1:nsites
        id = clusterid(uf, site)
        model.spins[site] = clusterspin[id]
        clustersize[id] += 1
    end
    return SWInfo(activated_bonds, clustersize, clusterspin)
end

function SW_update!(model::Potts, T::Real)
    p = -expm1(-1.0/T)
    nsites = numsites(model.lat)
    nbonds = numbonds(model.lat)
    activated_bonds = 0
    uf = UnionFind(nsites)
    @inbounds for bond in 1:nbonds
        s1,s2 = source(model.lat, bond), target(model.lat, bond)
        if model.spins[s1] == model.spins[s2] && rand() < p
            activated_bonds += 1
            unify!(uf, s1,s2)
        end
    end
    nc = clusterize!(uf)
    clustersize = ones(Int, nc)
    clusterspin = rand(1:model.Q, nc)

    @inbounds for site in 1:nsites
        id = clusterid(uf, site)
        model.spins[site] = clusterspin[id]
        clustersize[id] += 1
    end
    return SWInfo(activated_bonds, clustersize, clusterspin)
end

function SW_update!(model::XY, T::Real)
    nsites = numsites(model.lat)
    nbonds = numbonds(model.lat)
    m2b = -2/T
    m = 0.5*rand()
    pspins = zeros(nsites)
    @inbounds for s in 1:nsites
        pspins[s] = cospi(2(model.spins[s]-m))
    end
    uf = UnionFind(nsites)
    @inbounds for bond in 1:nbonds
        s1,s2 = source(model.lat, bond), target(model.lat, bond)
        if rand() < -expm1(m2b*pspins[s1]*pspins[s2])
            unify!(uf, s1,s2)
        end
    end
    nc = clusterize!(uf)
    toflips = rand([1,-1], nc)
    clustersize = zeros(Int, nc)
    @inbounds for site in 1:nsites
        id = clusterid(uf, site)
        s = model.spins[site]
        s = ifelse(toflips[id] > 0, s, 2m+0.5-s)
        s = mod(s, 1.0)
        clustersize[id] += 1
        model.spins[site] = s
    end
    return clustersize
end

