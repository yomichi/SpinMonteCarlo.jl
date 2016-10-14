"""
    measure(model::Ising, T::Real)
    measure(model::Potts, T::Real)

return magnetization density `M` and total energy `E`.
"""
function measure(model::Ising, T::Real)
    lat = model.lat
    nsites = numsites(lat)
    nbonds = numbonds(lat)

    M = mean(model.spins)
    E = 0.0
    @inbounds for b in 1:nbonds
        s1, s2 = source(lat, b), target(lat, b)
        E += ifelse(model.spins[s1] == model.spins[s2], -1.0, 1.0)
    end
    return M,E
end

function measure(model::Potts, T::Real)
    lat = model.lat
    nsites = numsites(lat)
    nbonds = numbonds(lat)
    invQ = 1.0/model.Q

    M = 0.0
    @inbounds for s in 1:nsites
        M += ifelse(model.spins[s]==1, 1.0-invQ, -invQ)
    end
    M /= nsites
    E = 0.0
    @inbounds for b in 1:nbonds
        s1, s2 = source(lat, b), target(lat, b)
        E -= ifelse(model.spins[s1] == model.spins[s2], 1.0, 0.0)
    end
    return M,E
end


"""
    measure(model::XY, T::Real)

return magnetization density `M`, total energy `E`, and helicity modulus `U`.
"""
function measure(model::XY, T::Real)
    lat = model.lat
    nsites = numsites(lat)
    nbonds = numbonds(lat)
    D = dim(lat)
    invN = 1.0/nsites
    beta = 1.0/T

    M = zeros(2)
    E = 0.0
    U1 = zeros(2)
    U2 = zeros(2)

    @inbounds for s in 1:nsites
        M[1] += cospi(2model.spins[s])
        M[2] += sinpi(2model.spins[s])
    end

    @inbounds for b in 1:nbonds
        i, j = source(lat, b), target(lat,b)
        dir = bonddirectory(lat, b)
        dt = mod(model.spins[j] - model.spins[i] + 2.0, 1.0)
        dt = ifelse(0.5 <= dt, dt - 1.0, dt)
        E -= cospi(2dt)

        for d in 1:D
            U1[d] += cospi(2dt) * dir[d]^2
            U2[d] += sinpi(2dt) * dir[d]
        end
    end
    for d in 1:D
        M[d] *= invN
        U1[d] -= beta * U2[d]^2
        U1[d] *= invN
    end
    return M, E, U1
end
