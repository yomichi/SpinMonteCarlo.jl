@doc raw"""
    postproc(model::AshkinTeller, param::Parameter, obs::MCObservableSet)

# Observables to be calculated
In the following, $m$ is total magnetization per site and $\epsilon$ is total energy per site.

- `"Binder Ratio"`
    - $R := \frac{\left \langle m^4 \right \rangle}{\left \langle m^2 \right\rangle^2}$
- `"Susceptibility"`
    - $\chi := \frac{N}{T}\left(\left\langle m^2\right\rangle\right)$
- `"Connected Susceptibility"`
    - $\frac{N}{T}\left(\left\langle m^2\right\rangle - \left\langle |m| \right\rangle^2\right)$
- `"Binder Ratio sigma"`
    - Binder ratio with respct to ``\sigma``
- `"Susceptibility sigma"`
    - Susceptibility with respct to ``\sigma``
- `"Connected Susceptibility sigma"`
    - Connected susceptibility with respct to ``\sigma``
- `"Binder Ratio tau"`
    - Binder ratio with respct to ``\tau``
- `"Susceptibility tau"`
    - Susceptibility with respct to ``\tau``
- `"Connected Susceptibility tau"`
    - Connected susceptibility with respct to ``\tau``
- `"Specific Heat"`
    - $\frac{N}{T^2}\left(\left\langle \epsilon^2\right\rangle - \left\langle \epsilon \right\rangle^2\right)$
"""
function postproc(model::AshkinTeller, param::Parameter, obs::MCObservableSet)
    nsites = numsites(model)
    T = convert(Float64, param["T"])
    beta = 1.0/T

    jk = jackknife(obs)

    jk["Binder Ratio"] = jk["|Magnetization|^4"] / (jk["|Magnetization|^2"]^2)
    jk["Susceptibility"] = (nsites*beta)*jk["|Magnetization|^2"]
    jk["Connected Susceptibility"] = (nsites*beta)*(jk["|Magnetization|^2"] - jk["|Magnetization|"]^2)

    jk["Binder Ratio sigma"] = jk["Magnetization sigma^4"] / (jk["Magnetization sigma^2"]^2)
    jk["Susceptibility sigma"] = (nsites*beta)*jk["Magnetization sigma^2"]
    jk["Connected Susceptibility sigma"] = (nsites*beta)*(jk["Magnetization sigma^2"] - jk["|Magnetization sigma|"]^2)

    jk["Binder Ratio tau"] = jk["Magnetization tau^4"] / (jk["Magnetization tau^2"]^2)
    jk["Susceptibility tau"] = (nsites*beta)*jk["Magnetization tau^2"]
    jk["Connected Susceptibility tau"] = (nsites*beta)*(jk["Magnetization tau^2"] - jk["|Magnetization tau|"]^2)

    jk["Specific Heat"] = (nsites*beta*beta)*(jk["Energy^2"] - jk["Energy"]^2)
    return jk
end
