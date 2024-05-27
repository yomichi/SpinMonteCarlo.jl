function simple_estimator(model::XY, T::Real, Js::AbstractArray, _=nothing)
    nsites = numsites(model)
    nbonds = numbonds(model)
    D = dim(model)
    invN = 1.0 / nsites
    beta = 1.0 / T

    M = zeros(2)
    E = 0.0
    U1 = zeros(2)
    U2 = zeros(2)

    @inbounds for s in 1:nsites
        M[1] += cospi(2model.spins[s])
        M[2] += sinpi(2model.spins[s])
    end

    @inbounds for b in bonds(model)
        i, j = source(b), target(b)
        dir = bonddirection(b)
        dt = mod(model.spins[j] - model.spins[i] + 2.0, 1.0)
        dt = ifelse(0.5 <= dt, dt - 1.0, dt)
        E -= cospi(2dt) * Js[bondtype(b)]

        for d in 1:D
            U1[d] += cospi(2dt) * dir[d]^2
            U2[d] += sinpi(2dt) * dir[d]
        end
    end
    for d in 1:2
        M[d] *= invN
        U1[d] -= beta * U2[d]^2
        U1[d] *= invN
    end
    E /= nsites

    res = Measurement()
    for (i, c) in enumerate(["x", "y"])
        res["Magnetization $c"] = M[i]
        res["|Magnetization $c|"] = abs(M[i])
        res["Magnetization $c^2"] = M[i]^2
        res["Magnetization $c^4"] = M[i]^4
        res["Helicity Modulus $c"] = U1[i]
    end
    M2 = sum(abs2, M)
    res["|Magnetization|"] = sqrt(M2)
    res["|Magnetization|^2"] = M2
    res["|Magnetization|^4"] = M2^2
    res["Energy"] = E
    res["Energy^2"] = E^2

    return res
end
