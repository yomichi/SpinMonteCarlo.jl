"""
    measure(model::XY, T::Real)

return magnetization `M`, energy `E`, and helicity modulus `U`.
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
