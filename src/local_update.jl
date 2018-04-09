"""
    local_update!(model, T::Real, J::Real; measure::Bool=true)
    local_update!(model, T::Real, Js::AbstractArray; measure::Bool=true)
update spin configuration by local spin flip and Metropolice algorithm under the temperature `T`
"""
@inline function local_update!(model::Model, T::Real, J::Real; measure::Bool=true)
    Js = J*ones(numbondtypes(model))
    return local_update!(model,T,Js,measure=measure)
end
function local_update!(model::Ising, T::Real, Js::AbstractArray; measure::Bool=true)
    nsites = numsites(model)
    nbonds = numbonds(model)
    mbeta = -1.0/T

    @inbounds for site in 1:nsites
        center = model.spins[site]
        de = 0.0
        for (n,b) in neighbors(model, site)
            de += 2center * model.spins[n] * Js[bondtype(model,b)]
        end
        if rand() < exp(mbeta*de)
            model.spins[site] *= -1
        end
    end

    res = Measurement()
    if measure
        M, E = simple_estimate(model, T, Js)
        res["M"] = M
        res["M2"] = M^2
        res["M4"] = M^4
        res["E"] = E
        res["E2"] = E^2
    end

    return res
end

function local_update!(model::Potts, T::Real, Js::AbstractArray; measure::Bool=true)
    nsites = numsites(model)
    nbonds = numbonds(model)
    mbeta = -1.0/T

    @inbounds for site in 1:nsites
        center = model.spins[site]
        new_center = mod1(center+rand(1:(model.Q-1)), model.Q)
        de = 0.0
        for (n,b) in neighbors(model, site)
            de += ifelse(center == model.spins[n], Js[bondtype(model, b)], 0.0)
            de -= ifelse(new_center == model.spins[n], Js[bondtype(model, b)], 0.0)
        end
        if rand() < exp(mbeta*de)
            model.spins[site] = new_center
        end
    end

    res = Measurement()
    if measure
        M, E = simple_estimate(model, T, Js)
        res["M"] = M
        res["M2"] = M^2
        res["M4"] = M^4
        res["E"] = E
        res["E2"] = E^2
    end

    return res
end

function local_update!(model::Clock, T::Real, Js::AbstractArray; measure::Bool=true)
    nsites = numsites(model)
    nbonds = numbonds(model)
    mbeta = -1.0/T
    iQ2 = 2.0/model.Q

    @inbounds for site in 1:nsites
        center = model.spins[site]
        new_center = mod1(center+rand(1:(model.Q-1)), model.Q)
        de = 0.0
        for (n,b) in neighbors(model, site)
            de += model.cosines[mod1(model.spins[n]-center, model.Q)] * Js[bondtype(model,b)]
            de -= model.cosines[mod1(model.spins[n]-new_center, model.Q)] * Js[bondtype(model,b)]
        end
        if rand() < exp(mbeta*de)
            model.spins[site] = new_center
        end
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

function local_update!(model::XY, T::Real, Js::AbstractArray; measure::Bool=true)
    nsites = numsites(model)
    nbonds = numbonds(model)
    mbeta = -1.0/T

    @inbounds for site in 1:nsites
        center = model.spins[site]
        new_center = rand()
        de = 0.0
        for (n,b) in neighbors(model, site)
            de += cospi(2(center - model.spins[n])) * Js[bondtype(model,b)]
            de -= cospi(2(new_center - model.spins[n])) * Js[bondtype(model,b)]
        end
        if rand() < exp(mbeta*de)
            model.spins[site] = new_center
        end
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

