@doc raw"""
    postproc(model::Potts, param::Parameter, obs::MCObservableSet)

# Observables to be calculated
In the following, $m$ is total magnetization per site and $\epsilon$ is total energy per site.

- `"Binder Ratio"`
    - $R := \frac{\left \langle m^4 \right \rangle}{\left \langle m^2 \right\rangle^2}$
- `"Susceptibility"`
    - $\chi := \frac{N}{T}\left(\left\langle m^2\right\rangle\right)$
- `"Connected Susceptibility"`
    - $\frac{N}{T}\left(\left\langle m^2\right\rangle - \left\langle |m| \right\rangle^2\right)$
- `"Specific Heat"`
    - $\frac{N}{T^2}\left(\left\langle \epsilon^2\right\rangle - \left\langle \epsilon \right\rangle^2\right)$
"""
function postproc(model::Potts, param::Parameter, obs::MCObservableSet)
    nsites = numsites(model)
    T = convert(Float64, param["T"])
    beta = 1.0/T

    jk = jackknife(obs)
    jk["Binder Ratio"] = jk["Magnetization^4"] / (jk["Magnetization^2"]^2)
    jk["Susceptibility"] = (nsites*beta)*jk["Magnetization^2"]
    jk["Connected Susceptibility"] = (nsites*beta)*(jk["Magnetization^2"] - jk["|Magnetization|"]^2)
    jk["Specific Heat"] = (nsites*beta*beta)*(jk["Energy^2"] - jk["Energy"]^2)
    return jk
end
