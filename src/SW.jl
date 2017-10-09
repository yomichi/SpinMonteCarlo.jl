type SWInfo
    activated_bonds :: Vector{Int}
    clustersize :: Vector{Int}
    clusterspin :: Vector{Int}
end

numclusters(sw::SWInfo) = length(sw.clustersize)

"""
    SW_update!(model, T::Real, J::Real; measure::Bool=true)
    SW_update!(model, T::Real, Js::AbstractArray; measure::Bool=true)
    
update spin configuration by Swendsen-Wang algorithm under the temperature `T`.
"""
function SW_update!(model::Model, T::Real, J::Real; measure::Bool=true)
    Js = J*ones(numbondtypes(model))
    return SW_update!(model, T, Js, measure=measure)
end
function SW_update!(model::Ising, T::Real, Js::AbstractArray; measure::Bool=true)
    ps = -expm1.((-2.0/T).*Js)
    nsites = numsites(model)
    nbonds = numbonds(model)
    nbt = numbondtypes(model)
    activated_bonds = zeros(Int,nbt)
    uf = UnionFind(nsites)
    @inbounds for bond in 1:nbonds
        s1,s2 = source(model, bond), target(model, bond)
        bt = bondtype(model,bond)
        if model.spins[s1] == model.spins[s2] && rand() < ps[bt]
            activated_bonds[bt] += 1
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
        M, M2, M4, E, E2 = improved_estimate(model, T, Js, swinfo)
        res["M"] = M
        res["M2"] = M2
        res["M4"] = M4
        res["E"] = E
        res["E2"] = E2
    end
    return res
end

function SW_update!(model::Potts, T::Real, Js::AbstractArray; measure::Bool=true)
    ps = -expm1.((-1.0/T).*Js)
    nsites = numsites(model)
    nbonds = numbonds(model)
    nbt = numbondtypes(model)
    activated_bonds = zeros(Int,nbt)
    uf = UnionFind(nsites)
    @inbounds for bond in 1:nbonds
        s1,s2 = source(model, bond), target(model, bond)
        bt = bondtype(model,bond)
        if model.spins[s1] == model.spins[s2] && rand() < ps[bt]
            activated_bonds[bt] += 1
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
        M, M2, M4, E, E2 = improved_estimate(model, T, Js, swinfo)
        res["M"] = M
        res["M2"] = M2
        res["M4"] = M4
        res["E"] = E
        res["E2"] = E2
    end
    return res
end

function SW_update!(model::Clock, T::Real, Js::AbstractArray; measure::Bool=true)
    nsites = numsites(model)
    nbonds = numbonds(model)
    m2bJ = (-2/T).*Js
    m = rand(1:model.Q)-1
    rspins = zeros(Int, nsites)
    @inbounds for s in 1:nsites
        rspins[s] = mod1(model.spins[s]-m, model.Q)
    end
    uf = UnionFind(nsites)
    @inbounds for bond in 1:nbonds
        s1,s2 = source(model, bond), target(model, bond)
        bt = bondtype(model,bond)
        if rand() < -expm1(m2bJ[bt]*model.sines_sw[rspins[s1]]*model.sines_sw[rspins[s2]])
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
        M, E, U = simple_estimate(model, T, Js)
        M2 = sum(abs2,M)
        res["M"] = M
        res["M2"] = M2
        res["M4"] = M2^2
        res["E"] = E
        res["E2"] = E^2
        res["U"] = U
    end

    return res
end

function SW_update!(model::XY, T::Real, Js::AbstractArray; measure::Bool=true)
    nsites = numsites(model)
    nbonds = numbonds(model)
    m2bJ = (-2/T).*Js
    m = 0.5*rand()
    pspins = zeros(nsites)
    @inbounds for s in 1:nsites
        pspins[s] = cospi(2(model.spins[s]-m))
    end
    uf = UnionFind(nsites)
    @inbounds for bond in 1:nbonds
        s1,s2 = source(model, bond), target(model, bond)
        bt = bondtype(model,bond)
        if rand() < -expm1(m2bJ[bt]*pspins[s1]*pspins[s2])
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
        M, E, U = simple_estimate(model, T, Js)
        M2 = sum(abs2,M)
        res["M"] = M
        res["M2"] = M2
        res["M4"] = M2^2
        res["E"] = E
        res["E2"] = E^2
        res["U"] = U
    end

    return res
end

