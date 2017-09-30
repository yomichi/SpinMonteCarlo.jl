type SWInfo
    activated_bonds :: Int
    clustersize :: Vector{Int}
    clusterspin :: Vector{Int}
end

numclusters(sw::SWInfo) = length(sw.clustersize)

"""
    SW_update!(model, T::Real; measure::Bool=true)
    
update spin configuration by Swendsen-Wang algorithm under the temperature `T`.
"""
function SW_update!(model::Ising, T::Real; measure::Bool=true)
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
    clustersize = zeros(Int, nc)
    clusterspin = rand([1,-1], nc)

    @inbounds for site in 1:nsites
        id = clusterid(uf, site)
        model.spins[site] = clusterspin[id]
        clustersize[id] += 1
    end
    swinfo = SWInfo(activated_bonds, clustersize, clusterspin)

    res = Measurement()
    if measure
        M, M2, M4, E, E2 = improved_estimate(model, T, swinfo)
        res[:M] = M
        res[:M2] = M2
        res[:M4] = M4
        res[:E] = E
        res[:E2] = E2
    end
    return res
end

function SW_update!(model::Potts, T::Real; measure::Bool=true)
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
    clustersize = zeros(Int, nc)
    clusterspin = rand(1:model.Q, nc)

    @inbounds for site in 1:nsites
        id = clusterid(uf, site)
        model.spins[site] = clusterspin[id]
        clustersize[id] += 1
    end
    swinfo =  SWInfo(activated_bonds, clustersize, clusterspin)

    res = Measurement()
    if measure
        M, M2, M4, E, E2 = improved_estimate(model, T, swinfo)
        res[:M] = M
        res[:M2] = M2
        res[:M4] = M4
        res[:E] = E
        res[:E2] = E2
    end
    return res
end

function SW_update!(model::Clock, T::Real; measure::Bool=true)
    nsites = numsites(model.lat)
    nbonds = numbonds(model.lat)
    m2b = -2/T
    m = rand(1:model.Q)
    rspins = zeros(Int, nsites)
    @inbounds for s in 1:nsites
        rspins[s] = mod1(model.spins[s]-m, model.Q)
    end
    uf = UnionFind(nsites)
    @inbounds for bond in 1:nbonds
        s1,s2 = source(model.lat, bond), target(model.lat, bond)
        if rand() < -expm1(m2b*model.sines_sw[rspins[s1]]*model.sines_sw[rspins[s2]])
            unify!(uf, s1,s2)
        end
    end
    nc = clusterize!(uf)
    toflips = rand(Bool, nc)
    clustersize = zeros(Int, nc)
    @inbounds for site in 1:nsites
        id = clusterid(uf, site)
        s = ifelse(toflips[id], model.Q-rspins[site]+1, rspins[site])
        clustersize[id] += 1
        model.spins[site] = mod1(s+m, model.Q)
    end

    res = Measurement()
    if measure
        M, E, U = simple_estimate(model, T)
        M2 = sum(abs2,M)
        res[:M] = M
        res[:M2] = M2
        res[:M4] = M2^2
        res[:E] = E
        res[:E2] = E^2
        res[:U] = U
    end

    return res
end

function SW_update!(model::XY, T::Real; measure::Bool=true)
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
    toflips = rand(Bool, nc)
    clustersize = zeros(Int, nc)
    @inbounds for site in 1:nsites
        id = clusterid(uf, site)
        s = model.spins[site]
        s = ifelse(toflips[id], s, 2m+0.5-s)
        s = mod(s, 1.0)
        clustersize[id] += 1
        model.spins[site] = s
    end

    res = Measurement()
    if measure
        M, E, U = simple_estimate(model, T)
        M2 = sum(abs2,M)
        res[:M] = M
        res[:M2] = M2
        res[:M4] = M2^2
        res[:E] = E
        res[:E2] = E^2
        res[:U] = U
    end

    return res
end

