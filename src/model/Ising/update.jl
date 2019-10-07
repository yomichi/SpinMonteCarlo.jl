function local_update!(model::Ising, T::Real, Js::AbstractArray)
    rng = model.rng
    nsites = numsites(model)
    nbonds = numbonds(model)
    mbeta = -1.0/T

    @inbounds for site in 1:nsites
        center = model.spins[site]
        de = 0.0
        for (n,b) in neighbors(model, site)
            de += 2center * model.spins[n] * Js[bondtype(model,b)]
        end
        if rand(rng) < exp(mbeta*de)
            model.spins[site] *= -1
        end
    end

    return nothing
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

function Wolff_update!(model::Ising, T::Real, Js::AbstractArray)
    rng = model.rng
    ps = -expm1.((-2.0/T).*Js)
    nsites = numsites(model)

    clustersize = 0
    st = Stack(Deque{Int}())
    center = rand(rng, 1:nsites)
    sp = model.spins[center]
    model.spins[center] *= -1
    push!(st, center)
    @inbounds while !isempty(st)
        clustersize += 1
        s = pop!(st)
        for (n,b) in neighbors(model, s)
            bt = bondtype(model,b)
            if model.spins[n] == sp && rand(rng) < ps[bt]
                model.spins[n] *= -1
                push!(st, n)
            end
        end
    end

    return nothing
end

@gen_convert_parameter(Ising, ("T", 1, 1.0), ("J", numbondtypes, 1.0))
