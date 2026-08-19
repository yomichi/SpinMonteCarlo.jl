function local_update!(model::Ising, T::Real, Js::AbstractArray)
    rng = model.rng
    nsites = numsites(model)
    nbonds = numbonds(model)
    mbeta = -1.0 / T

    @inbounds for _ in 1:nsites
        site = rand(rng, 1:nsites)
        center = model.spins[site]
        de = 0.0
        for (n, b) in neighbors(model, site)
            de += 2center * model.spins[n] * Js[bondtype(model, b)]
        end
        if rand(rng) < exp(mbeta * de)
            model.spins[site] *= -1
        end
    end

    return nothing
end

function SW_update!(model::Ising, T::Real, Js::AbstractArray)
    rng = model.rng
    ps = -expm1.((-2.0 / T) .* abs.(Js))
    nsites = numsites(model)
    nbonds = numbonds(model)
    nbt = numbondtypes(model)
    activated_bonds = zeros(Int, nbt)
    uf = UnionFind(nsites)
    @inbounds for bond in bonds(model)
        s1, s2 = source(bond), target(bond)
        bt = bondtype(bond)
        if Js[bt] * model.spins[s1] * model.spins[s2] > 0 && rand(rng) < ps[bt]
            activated_bonds[bt] += 1
            unify!(uf, s1, s2)
        end
    end
    nc = clusterize!(uf)
    clustersize = zeros(Int, nc)
    clustermag = zeros(Int, nc)
    clusterspin = rand(rng, [1, -1], nc)

    @inbounds for site in 1:nsites
        id = clusterid(uf, site)
        clustersize[id] += 1
        clustermag[id] += model.spins[site]
        model.spins[site] *= clusterspin[id]
    end
    return SWInfo(activated_bonds, clustersize, clusterspin, clustermag)
end

function Wolff_update!(model::Ising, T::Real, Js::AbstractArray)
    rng = model.rng
    nsites = numsites(model)

    in_cluster = falses(nsites)
    cluster = Int[]
    st = Int[]
    center = rand(rng, 1:nsites)
    in_cluster[center] = true
    push!(cluster, center)
    push!(st, center)
    @inbounds while !isempty(st)
        s = pop!(st)
        for (n, b) in neighbors(model, s)
            if in_cluster[n]
                continue
            end
            bt = bondtype(model, b)
            x = (-2.0 / T) * Js[bt] * model.spins[s] * model.spins[n]
            p = -expm1(min(0.0, x))
            if rand(rng) < p
                in_cluster[n] = true
                push!(cluster, n)
                push!(st, n)
            end
        end
    end

    @inbounds for site in cluster
        model.spins[site] *= -1
    end

    return nothing
end

@gen_convert_parameter(Ising, ("T", 1, 1.0), ("J", numbondtypes, 1.0))
