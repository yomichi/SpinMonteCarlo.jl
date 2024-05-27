function local_update!(model::Clock, T::Real, Js::AbstractArray)
    rng = model.rng
    nsites = numsites(model)
    nbonds = numbonds(model)
    mbeta = -1.0 / T
    iQ2 = 2.0 / model.Q

    @inbounds for site in 1:nsites
        center = model.spins[site]
        new_center = mod1(center + rand(rng, 1:(model.Q - 1)), model.Q)
        de = 0.0
        for (n, b) in neighbors(model, site)
            de += model.cosines[mod1(model.spins[n] - center, model.Q)] *
                  Js[bondtype(model, b)]
            de -= model.cosines[mod1(model.spins[n] - new_center, model.Q)] *
                  Js[bondtype(model, b)]
        end
        if rand(rng) < exp(mbeta * de)
            model.spins[site] = new_center
        end
    end

    return nothing
end

function SW_update!(model::Clock, T::Real, Js::AbstractArray)
    rng = model.rng
    nsites = numsites(model)
    nbonds = numbonds(model)
    m2bJ = (-2 / T) .* Js
    m = rand(rng, 1:(model.Q)) - 1
    rspins = zeros(Int, nsites)
    @inbounds for s in 1:nsites
        rspins[s] = mod1(model.spins[s] - m, model.Q)
    end
    uf = UnionFind(nsites)
    @inbounds for bond in bonds(model)
        s1, s2 = source(bond), target(bond)
        bt = bondtype(bond)
        if rand(rng) <
           -expm1(m2bJ[bt] * model.sines_sw[rspins[s1]] * model.sines_sw[rspins[s2]])
            unify!(uf, s1, s2)
        end
    end
    nc = clusterize!(uf)
    toflips = rand(rng, Bool, nc)
    clustersize = zeros(Int, nc)
    @inbounds for site in 1:nsites
        id = clusterid(uf, site)
        s = ifelse(toflips[id], model.Q - rspins[site] + 1, rspins[site])
        clustersize[id] += 1
        model.spins[site] = mod1(s + m, model.Q)
    end

    return nothing
end

function Wolff_update!(model::Clock, T::Real, Js::AbstractArray)
    rng = model.rng
    nsites = numsites(model)
    b2J = (2.0 / T) .* Js

    clustersize = 0
    st = Stack(Deque{Int}())

    m = rand(rng, 1:(model.Q)) - 1
    center = rand(rng, 1:nsites)
    push!(st, center)
    r = mod1(model.spins[center] - m, model.Q)
    model.spins[center] = mod1(m - r + 1, model.Q)

    @inbounds while !isempty(st)
        clustersize += 1
        center = pop!(st)
        cr = mod1(model.spins[center] - m, model.Q)
        for (n, b) in neighbors(model, center)
            nr = mod1(model.spins[n] - m, model.Q)
            bt = bondtype(model, b)
            ## Note: center already flipped
            if rand(rng) < -expm1(b2J[bt] * model.sines_sw[cr] * model.sines_sw[nr])
                push!(st, n)
                model.spins[n] = mod1(m - nr + 1, model.Q)
            end
        end
    end

    return nothing
end

@gen_convert_parameter(Clock, ("T", 1, 1.0), ("J", numbondtypes, 1.0))
