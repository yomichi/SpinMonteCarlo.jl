function local_update!(model::XY, T::Real, Js::AbstractArray)
    rng = model.rng
    nsites = numsites(model)
    nbonds = numbonds(model)
    mbeta = -1.0/T

    @inbounds for site in 1:nsites
        center = model.spins[site]
        new_center = rand(rng)
        de = 0.0
        for (n,b) in neighbors(model, site)
            de += cospi(2(center - model.spins[n])) * Js[bondtype(model,b)]
            de -= cospi(2(new_center - model.spins[n])) * Js[bondtype(model,b)]
        end
        if rand(rng) < exp(mbeta*de)
            model.spins[site] = new_center
        end
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

function Wolff_update!(model::XY, T::Real, Js::AbstractArray)
    rng = model.rng
    nsites = numsites(model)
    b2J = (2.0/T).*Js

    clustersize = 0
    st = Stack(Deque{Int}())
    
    m = 0.5*rand(rng)
    center = rand(rng,1:nsites)
    push!(st, center)

    model.spins[center] = mod(2m+0.5-model.spins[center], 1.0)

    @inbounds while !isempty(st)
        clustersize += 1
        center = pop!(st)
        cp = cospi(2(model.spins[center]-m))

        for (n,b) in neighbors(model, center)
            np = cospi(2(model.spins[n]-m))
            bt = bondtype(model,b)
            ## Note: center already flipped
            if rand(rng) < -expm1(b2J[bt]*cp*np)
                push!(st, n)
                model.spins[n] = mod(2m+0.5-model.spins[n], 1.0)
            end
        end
    end

    return nothing
end

@gen_convert_parameter(XY, ("T", 1, 1.0), ("J", numbondtypes, 1.0))
