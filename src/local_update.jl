@doc """
    local_update!(model, param)
    local_update!(model, T::Real, Js::AbstractArray)

Updates spin configuration by local spin flip and Metropolice algorithm 
under the temperature `T = param["T"]` and coupling constants `J = param["J"]`
"""
@inline function local_update!(model::Model, param::Parameter)
    p = convert_parameter(model, param)
    return local_update!(model, p...)
end

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

function local_update!(model::Potts, T::Real, Js::AbstractArray)
    rng = model.rng
    nsites = numsites(model)
    nbonds = numbonds(model)
    mbeta = -1.0/T

    @inbounds for site in 1:nsites
        center = model.spins[site]
        new_center = mod1(center+rand(rng, 1:(model.Q-1)), model.Q)
        de = 0.0
        for (n,b) in neighbors(model, site)
            de += ifelse(center == model.spins[n], Js[bondtype(model, b)], 0.0)
            de -= ifelse(new_center == model.spins[n], Js[bondtype(model, b)], 0.0)
        end
        if rand(rng) < exp(mbeta*de)
            model.spins[site] = new_center
        end
    end

    return nothing
end

function local_update!(model::Clock, T::Real, Js::AbstractArray)
    rng = model.rng
    nsites = numsites(model)
    nbonds = numbonds(model)
    mbeta = -1.0/T
    iQ2 = 2.0/model.Q

    @inbounds for site in 1:nsites
        center = model.spins[site]
        new_center = mod1(center+rand(rng, 1:(model.Q-1)), model.Q)
        de = 0.0
        for (n,b) in neighbors(model, site)
            de += model.cosines[mod1(model.spins[n]-center, model.Q)] * Js[bondtype(model,b)]
            de -= model.cosines[mod1(model.spins[n]-new_center, model.Q)] * Js[bondtype(model,b)]
        end
        if rand(rng) < exp(mbeta*de)
            model.spins[site] = new_center
        end
    end

    return nothing
end

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

