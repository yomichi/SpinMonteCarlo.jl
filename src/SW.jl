@doc doc"""
Information of clusters in Swendsen-Wang algorithm.

# Fields
- `activated_bonds` : The number of activated (connected) bonds of each cluster.
- `clustersize` : The number of sites in each cluster.
- `clusterspin` : Spin variable of each cluster (e.g., 1 or -1 for `Ising`).
"""
mutable struct SWInfo
    activated_bonds :: Vector{Int}
    clustersize :: Vector{Int}
    clusterspin :: Vector{Int}
end

numclusters(sw::SWInfo) = length(sw.clustersize)

@doc """
    SW_update!(model, param::Parameter)
    SW_update!(model, T::Real, Js::AbstractArray)
    
Updates spin configuration by Swendsen-Wang algorithm
under temperature `T=param["T"]` and coupling constants `J=param["J"]`
"""
@inline function SW_update!(model::Union{Ising, Potts, Clock, XY}, param::Parameter)
    p = convert_parameter(model, param)
    return SW_update!(model, p...)
end

function SW_update!(model::Ising, T::Real, Js::AbstractArray)
    rng = model.rng
    ps = -expm1.((-2.0/T).*Js)
    nsites = numsites(model)
    nbonds = numbonds(model)
    nbt = numbondtypes(model)
    activated_bonds = zeros(Int,nbt)
    uf = UnionFind(nsites)
    @inbounds for bond in bonds(model)
        s1,s2 = source(bond), target(bond)
        bt = bondtype(bond)
        if model.spins[s1] == model.spins[s2] && rand(rng) < ps[bt]
            activated_bonds[bt] += 1
            unify!(uf, s1,s2)
        end
    end
    nc = clusterize!(uf)
    clustersize = zeros(Int, nc)
    clusterspin = rand(rng,[1,-1], nc)

    @inbounds for site in 1:nsites
        id = clusterid(uf, site)
        model.spins[site] = clusterspin[id]
        clustersize[id] += 1
    end
    return SWInfo(activated_bonds, clustersize, clusterspin)
end

function SW_update!(model::Potts, T::Real, Js::AbstractArray)
    rng = model.rng
    ps = -expm1.((-1.0/T).*Js)
    nsites = numsites(model)
    nbonds = numbonds(model)
    nbt = numbondtypes(model)
    activated_bonds = zeros(Int,nbt)
    uf = UnionFind(nsites)
    @inbounds for bond in bonds(model)
        s1,s2 = source(bond), target(bond)
        bt = bondtype(bond)
        if model.spins[s1] == model.spins[s2] && rand(rng) < ps[bt]
            activated_bonds[bt] += 1
            unify!(uf, s1,s2)
        end
    end
    nc = clusterize!(uf)
    clustersize = zeros(Int, nc)
    clusterspin = rand(rng,1:model.Q, nc)

    @inbounds for site in 1:nsites
        id = clusterid(uf, site)
        model.spins[site] = clusterspin[id]
        clustersize[id] += 1
    end
    return SWInfo(activated_bonds, clustersize, clusterspin)
end

function SW_update!(model::Clock, T::Real, Js::AbstractArray)
    rng = model.rng
    nsites = numsites(model)
    nbonds = numbonds(model)
    m2bJ = (-2/T).*Js
    m = rand(rng, 1:model.Q)-1
    rspins = zeros(Int, nsites)
    @inbounds for s in 1:nsites
        rspins[s] = mod1(model.spins[s]-m, model.Q)
    end
    uf = UnionFind(nsites)
    @inbounds for bond in bonds(model)
        s1,s2 = source(bond), target(bond)
        bt = bondtype(bond)
        if rand(rng) < -expm1(m2bJ[bt]*model.sines_sw[rspins[s1]]*model.sines_sw[rspins[s2]])
            unify!(uf, s1,s2)
        end
    end
    nc = clusterize!(uf)
    toflips = rand(rng, Bool, nc)
    clustersize = zeros(Int, nc)
    @inbounds for site in 1:nsites
        id = clusterid(uf, site)
        s = ifelse(toflips[id], model.Q-rspins[site]+1, rspins[site])
        clustersize[id] += 1
        model.spins[site] = mod1(s+m, model.Q)
    end

    return nothing
end

function SW_update!(model::XY, T::Real, Js::AbstractArray)
    rng = model.rng
    nsites = numsites(model)
    nbonds = numbonds(model)
    m2bJ = (-2/T).*Js
    m = 0.5*rand(rng)
    pspins = zeros(nsites)
    @inbounds for s in 1:nsites
        pspins[s] = cospi(2(model.spins[s]-m))
    end
    uf = UnionFind(nsites)
    @inbounds for bond in bonds(model)
        s1,s2 = source(bond), target(bond)
        bt = bondtype(bond)
        if rand(rng) < -expm1(m2bJ[bt]*pspins[s1]*pspins[s2])
            unify!(uf, s1,s2)
        end
    end
    nc = clusterize!(uf)
    toflips = rand(rng, Bool, nc)
    clustersize = zeros(Int, nc)
    @inbounds for site in 1:nsites
        id = clusterid(uf, site)
        s = model.spins[site]
        s = ifelse(toflips[id], s, 2m+0.5-s)
        s = mod(s, 1.0)
        clustersize[id] += 1
        model.spins[site] = s
    end

    return nothing
end

