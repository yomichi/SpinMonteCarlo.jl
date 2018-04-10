"""
    Wolff_update!(model, T::Real, J::Real)
    Wolff_update!(model, T::Real, Js::AbstractArray)

update spin configuration by Wolff algorithm under the temperature `T`.
"""
@inline function Wolff_update!(model::Model, T::Real, J::Real; measure::Bool=true)
    Js = J * ones(numbondtypes(model))
    return Wolff_update!(model, T, Js, measure=measure)
end
function Wolff_update!(model::Ising, T::Real, Js::AbstractArray; measure::Bool=true)
    rng = model.rng
    ps = -expm1.((-2.0/T).*Js)
    nsites = numsites(model)

    clustersize = 0
    st = Stack(Int)
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

    res = Measurement()
    if measure
        M, E = simple_estimate(model, T, Js)
        res["M"] = M
        res["M2"] = clustersize/nsites
        res["M4"] = M^4
        res["E"] = E
        res["E2"] = E^2
    end

    return res

end

function Wolff_update!(model::Potts, T::Real, Js::AbstractArray; measure::Bool=true)
    rng = model.rng
    ps = -expm1.((-1.0/T).*Js)
    nsites = numsites(model)

    clustersize = 0
    st = Stack(Int)
    center = rand(rng, 1:nsites)
    sp = model.spins[center]
    newsp = mod1(sp+rand(rng, 1:(model.Q-1)), model.Q)
    model.spins[center] = newsp
    push!(st, center)
    @inbounds while !isempty(st)
        clustersize += 1
        s = pop!(st)
        for (n,b) in neighbors(model, s)
            bt = bondtype(model,b)
            if model.spins[n] == sp && rand(rng) < ps[bt]
                model.spins[n] = newsp
                push!(st, n)
            end
        end
    end

    res = Measurement()
    if measure
        I2 = (model.Q-1)/(model.Q*model.Q)
        M, E = simple_estimate(model, T, Js)
        res["M"] = M
        res["M2"] = clustersize*I2/nsites
        res["M4"] = M^4
        res["E"] = E
        res["E2"] = E^2
    end

    return res
end

function Wolff_update!(model::Clock, T::Real, Js::AbstractArray; measure::Bool=true)
    rng = model.rng
    nsites = numsites(model)
    b2J = (2.0/T).*Js

    clustersize = 0
    st = Stack(Int)
    
    m = rand(rng, 1:model.Q)-1
    center = rand(rng, 1:nsites)
    push!(st, center)
    r = mod1(model.spins[center]-m, model.Q)
    model.spins[center] = mod1(m-r+1, model.Q)

    @inbounds while !isempty(st)
        clustersize += 1
        center = pop!(st)
        cr = mod1(model.spins[center]-m, model.Q)
        for (n,b) in neighbors(model, center)
            nr = mod1(model.spins[n]-m, model.Q)
            bt = bondtype(model,b)
            ## Note: center already flipped
            if rand(rng) < -expm1(b2J[bt]*model.sines_sw[cr]*model.sines_sw[nr])
                push!(st, n)
                model.spins[n] = mod1(m-nr+1,model.Q)
            end
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

function Wolff_update!(model::XY, T::Real, Js::AbstractArray; measure::Bool=true)
    rng = model.rng
    nsites = numsites(model)
    b2J = (2.0/T).*Js

    clustersize = 0
    st = Stack(Int)
    
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

