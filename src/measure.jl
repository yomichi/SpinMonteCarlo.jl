function measure(model::XY)
    lat = model.lat
    nsites = numsites(lat)
    nbonds = numbonds(lat)

    M = zeros(2)
    E = 0.0
    winding = zeros(2)

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

    end

end
