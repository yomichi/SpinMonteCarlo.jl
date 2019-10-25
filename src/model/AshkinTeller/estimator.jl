@doc raw"""
    simple_estimator(model::AshkinTeller, T::Real, Jsigma, Jtau, K)

Returns the following observables as `Dict{String, Any}`

# Observables
- `"Energy"`
    - energy density
- `"Energy^2"`
    - square of energy density
- `"|Magnetization|"`
    - absolute value of magnetization density, ``|m| = \sqrt{ (\sum_i \sigma_i )^2 + (\sum_i \tau_i)^2 } / N``
- `"|Magnetization|^2"`
    - square of magnetization density
- `"|Magnetization|^4"`
    - quadruple of magnetization density
- `"Magnetization sigma"`
    - magnetization density (sigma spin)
- `"|Magnetization sigma|"`
- `"Magnetization sigma^2"`
- `"Magnetization sigma^4"`
- `"Magnetization tau"`
    - magnetization density (tau spin)
- `"|Magnetization tau|"`
- `"Magnetization tau^2"`
- `"Magnetization tau^4"`

"""
function simple_estimator(model::AshkinTeller, T::Real, Jsigma, Jtau, K, _=nothing)
    nsites = numsites(model)
    nbonds = numbonds(model)

    Ms = mean(model.spins, dims=2)
    E = 0.0
    @inbounds for b in bonds(model)
        s1, s2 = source(b), target(b)
        bt = bondtype(b)
        E += model.spins[1,s1] * model.spins[1,s2] * Jsigma[bt]
        E += model.spins[2,s1] * model.spins[2,s2] * Jtau[bt]
        E += model.spins[1,s1] * model.spins[2,s1] * model.spins[1,s2] * model.spins[2,s2] * K[bt]
    end
    E /= nsites

    res = Measurement()
    res["|Magnetization|^2"] = Ms[1]^2 + Ms[2]^2
    res["|Magnetization|^4"] = res["|Magnetization|^2"]^2
    res["|Magnetization|"] = sqrt(res["|Magnetization|^2"])

    res["Magnetization sigma"] = Ms[1]
    res["|Magnetization sigma|"] = abs(Ms[1])
    res["Magnetization sigma^2"] = Ms[1]^2
    res["Magnetization sigma^4"] = Ms[1]^4

    res["Magnetization tau"] = Ms[2]
    res["|Magnetization tau|"] = abs(Ms[2])
    res["Magnetization tau^2"] = Ms[2]^2
    res["Magnetization tau^4"] = Ms[2]^4

    res["Energy"] = E
    res["Energy^2"] = E^2
    return res
end

default_estimator(model::AshkinTeller, update) = simple_estimator
