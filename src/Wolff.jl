"""
    Wolff_update!(model, T::Real)

update spin configuration by Wolff algorithm under the temperature `T`.
"""
function Wolff_update!(model::Ising, T::Real; measure::Bool=true)
    p = -expm1(-2.0/T)
    nsites = numsites(model)

    clustersize = 0
    st = Stack(Int)
    center = rand(1:nsites)
    sp = model.spins[center]
    model.spins[center] *= -1
    push!(st, center)
    @inbounds while !isempty(st)
        clustersize += 1
        s = pop!(st)
        for (n,_) in neighbors(model, s)
            if model.spins[n] == sp && rand() < p
                model.spins[n] *= -1
                push!(st, n)
            end
        end
    end

    res = Measurement()
    if measure
        M, E = simple_estimate(model, T)
        res[:M] = M
        res[:M2] = clustersize/nsites
        res[:M4] = M^4
        res[:E] = E
        res[:E2] = E^2
    end

    return res

end

function Wolff_update!(model::Potts, T::Real; measure::Bool=true)
    p = -expm1(-1.0/T)
    nsites = numsites(model)

    clustersize = 0
    st = Stack(Int)
    center = rand(1:nsites)
    sp = model.spins[center]
    newsp = mod1(sp+rand(1:(model.Q-1)), model.Q)
    model.spins[center] = newsp
    push!(st, center)
    @inbounds while !isempty(st)
        clustersize += 1
        s = pop!(st)
        for (n,_) in neighbors(model, s)
            if model.spins[n] == sp && rand() < p
                model.spins[n] = newsp
                push!(st, n)
            end
        end
    end

    res = Measurement()
    if measure
        I2 = (model.Q-1)/(model.Q*model.Q)
        M, E = simple_estimate(model, T)
        res[:M] = M
        res[:M2] = clustersize*I2/nsites
        res[:M4] = M^4
        res[:E] = E
        res[:E2] = E^2
    end

    return res
end

function Wolff_update!(model::Clock, T::Real; measure::Bool=true)
    nsites = numsites(model)
    b2 = 2.0/T

    clustersize = 0
    st = Stack(Int)
    
    m = rand(1:model.Q)-1
    center = rand(1:nsites)
    push!(st, center)
    r = mod1(model.spins[center]-m, model.Q)
    model.spins[center] = mod1(m-r+1, model.Q)

    @inbounds while !isempty(st)
        clustersize += 1
        center = pop!(st)
        cr = mod1(model.spins[center]-m, model.Q)
        for (n,_) in neighbors(model, center)
            nr = mod1(model.spins[n]-m, model.Q)
            ## Note: center already flipped
            if rand() < -expm1(b2*model.sines_sw[cr]*model.sines_sw[nr])
                push!(st, n)
                model.spins[n] = mod1(m-nr+1,model.Q)
            end
        end
    end

    res = Measurement()
    if measure
        M, E, U = simple_estimate(model, T)
        M2 = sum(abs2,M)
        res[:M] = M
        res[:M2] = M2
        res[:M4] = M2^2
        res[:E] = E
        res[:E2] = E^2
        res[:U] = U
    end

    return res
end

function Wolff_update!(model::XY, T::Real; measure::Bool=true)
    nsites = numsites(model)
    beta2 = 2.0/T

    clustersize = 0
    st = Stack(Int)
    
    m = 0.5*rand()
    center = rand(1:nsites)
    push!(st, center)

    model.spins[center] = mod(2m+0.5-model.spins[center], 1.0)

    @inbounds while !isempty(st)
        clustersize += 1
        center = pop!(st)
        cp = cospi(2(model.spins[center]-m))

        for (n,_) in neighbors(model, center)
            np = cospi(2(model.spins[n]-m))
            ## Note: center already flipped
            if rand() < -expm1(beta2*cp*np)
                push!(st, n)
                model.spins[n] = mod(2m+0.5-model.spins[n], 1.0)
            end
        end
    end

    res = Measurement()
    if measure
        M, E, U = simple_estimate(model, T)
        M2 = sum(abs2,M)
        res[:M] = M
        res[:M2] = M2
        res[:M4] = M2^2
        res[:E] = E
        res[:E2] = E^2
        res[:U] = U
    end

    return res
end

