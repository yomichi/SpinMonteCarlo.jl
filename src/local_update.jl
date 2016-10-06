"""
    local_update!(model, T::Real)
update spin configuration by local spin flip and Metropolice algorithm under the temperature `T`, and return the difference of total magnetization and energy.
"""
function local_update!(model::Ising, T::Real)
    nsites = numsites(model.lat)
    nbonds = numbonds(model.lat)
    mbeta = -1.0/T

    DM = 0
    DE = 0.0
    @inbounds for site in 1:nsites
        center = model.spins[site]
        de = 0.0
        for n in neighbors(model.lat, site)
            de += 2center * model.spins[n]
        end
        if rand() < exp(mbeta*de)
            model.spins[site] *= -1.0
            DE += de
            DM += 2model.spins[site]
        end
    end
    return DM, DE
end

function local_update!(model::XY, T::Real)
    nsites = numsites(model.lat)
    nbonds = numbonds(model.lat)
    mbeta = 1.0/T

    DM = zeros(2)
    DE = 0.0
    @inbounds for site in 1:nsites
        center = model.spins[site]
        new_center = rand()
        de = 0.0
        for n in neighbors(model.lat, site)
            de += cospi(2(center - model.spins[n]))
            de -= cospi(2(new_center - model.spins[n]))
        end
        if rand() < exp(mbeta*de)
            DM -= [cospi(2center), sinpi(2center)]
            DM += [cospi(2new_center), sinpi(2new_center)]
            model.spins[site] = new_center
            DE += de
        end
    end
    return DM, DE
end

