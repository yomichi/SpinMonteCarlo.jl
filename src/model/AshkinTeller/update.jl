function local_update!(model::AshkinTeller, T::Real, Jsigma, Jtau, K)
    rng = model.rng
    nsites = numsites(model)
    nbonds = numbonds(model)
    mbeta = -1.0 / T

    @inbounds for kind in 1:2
        Js = ifelse(kind == 1, Jsigma, Jtau)
        dual_kind = 3 - kind
        for site in 1:nsites
            center = model.spins[kind, site]
            dual_center = model.spins[dual_kind, site]
            de = 0.0
            for (n, b) in neighbors(model, site)
                de += 2center * model.spins[kind, n] *
                      (Js[bondtype(model, b)]
                       +
                       K[bondtype(model, b)] * model.spins[dual_kind, n] * dual_center)
            end
            if rand(rng) < exp(mbeta * de)
                model.spins[kind, site] *= -1
            end
        end
    end

    return nothing
end

function SW_update!(model::AshkinTeller, T::Real, Jsigma, Jtau, K)
    rng = model.rng
    m2b = -2.0 / T
    spins = model.spins
    nsites = numsites(model)
    nbonds = numbonds(model)
    nbt = numbondtypes(model)
    ps = zeros(nbt, 2)
    @inbounds for kind in 1:2
        uf = UnionFind(nsites)
        Js = ifelse(kind == 1, Jsigma, Jtau)
        ps[:, 1] .= -expm1.(m2b .* (Js .+ K))
        ps[:, 2] .= -expm1.(m2b .* (Js .- K))
        dual_kind = 3 - kind
        for bond in bonds(model)
            s1, s2 = source(bond), target(bond)
            bt = bondtype(bond)
            if spins[kind, s1] == spins[kind, s2] &&
               rand(rng) <
               ps[bt, ifelse(spins[dual_kind, s1] == spins[dual_kind, s2], 1, 2)]
                unify!(uf, s1, s2)
            end
        end
        nc = clusterize!(uf)
        clusterspin = rand(rng, [1, -1], nc)

        @inbounds for site in 1:nsites
            id = clusterid(uf, site)
            model.spins[kind, site] = clusterspin[id]
        end
    end
    return nothing
end

function Wolff_update!(model::AshkinTeller, T::Real, Jsigma, Jtau, K)
    rng = model.rng
    m2beta = -2.0 / T
    nsites = numsites(model)
    ps = zeros(numbondtypes(model), 2)

    @inbounds for kind in 1:2
        Js = ifelse(kind == 1, Jsigma, Jtau)
        ps[:, 1] .= -expm1.(m2beta .* (Js .+ K))
        ps[:, 2] .= -expm1.(m2beta .* (Js .- K))
        dual_kind = 3 - kind
        clustersize = 0
        st = Stack(Deque{Int}())
        center = rand(rng, 1:nsites)
        sp = model.spins[kind, center]
        model.spins[kind, center] *= -1
        push!(st, center)
        while !isempty(st)
            clustersize += 1
            s = pop!(st)
            for (n, b) in neighbors(model, s)
                bt = bondtype(model, b)
                if (model.spins[kind, n] == sp &&
                    rand(rng) < ps[bt,
                                   ifelse(model.spins[dual_kind, n] == model.spins[dual_kind, s], 1, 2)])
                    model.spins[kind, n] *= -1
                    push!(st, n)
                end
            end
        end
    end

    return nothing
end

@gen_convert_parameter(AshkinTeller, ("T", 1, 1.0), ("Jsigma", numbondtypes, 1.0),
                       ("Jtau", numbondtypes, 1.0), ("K", numbondtypes, 0.0))
