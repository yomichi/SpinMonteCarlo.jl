@doc """
    Wolff_update!(model, param::Parameter)
    Wolff_update!(model, T::Real, Js::AbstractArray)

Updates spin configuration by Wolff algorithm
under temperature `T=param["T"]` and coupling constants `J=param["J"]`
"""
@inline function Wolff_update!(model::Union{Ising, Potts, Clock, XY}, param::Parameter)
    p = convert_parameter(model, param)
    return Wolff_update!(model, p...)
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

function Wolff_update!(model::Potts, T::Real, Js::AbstractArray)
    rng = model.rng
    ps = -expm1.((-1.0/T).*Js)
    nsites = numsites(model)

    clustersize = 0
    st = Stack(Deque{Int}())
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
    return nothing
end

function Wolff_update!(model::Clock, T::Real, Js::AbstractArray)
    rng = model.rng
    nsites = numsites(model)
    b2J = (2.0/T).*Js

    clustersize = 0
    st = Stack(Deque{Int}())
    
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

