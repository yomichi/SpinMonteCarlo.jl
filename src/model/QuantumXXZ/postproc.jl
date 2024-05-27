@doc raw"""
    postproc(model::QuantumXXZ, param::Parameter, obs::MCObservableSet)

# Observables to be calculated
In the following, $s$ is sign of weight, $m$ is total magnetization per site,
and $\epsilon$ is total energy per site.

- `"Magnetization"`
    - $\left\langle m s\right\rangle\Big/\left\langle s \right\rangle$
- `"|Magnetization|"`
    - $\left\langle |m| s\right\rangle\Big/\left\langle s \right\rangle$
- `"Magnetization^2"`
    - $\left\langle m^2 s\right\rangle\Big/\left\langle s \right\rangle$
- `"Magnetization^4"`
    - $\left\langle m^4 s\right\rangle\Big/\left\langle s \right\rangle$
- `"Energy"`
    - $\left\langle \epsilon s\right\rangle\Big/\left\langle s \right\rangle$
- `"Energy^2"`
    - $\left\langle \epsilon^2 s\right\rangle\Big/\left\langle s \right\rangle$
- `"Binder Ratio"`
    - $\frac{\left \langle m^4 \right \rangle}{\left \langle m^2 \right\rangle^2}$
- `"Susceptibility"`
    - $\frac{N}{T}\left(\left\langle m^2\right\rangle\right)$
- `"Connected Susceptibility"`
    - $\frac{N}{T}\left(\left\langle m^2\right\rangle - \left\langle |m| \right\rangle^2\right)$
- `"Specific Heat"`
    - $\frac{N}{T^2}\left(\left\langle \epsilon^2\right\rangle - \left\langle \epsilon \right\rangle^2\right)$
"""
function postproc(model::QuantumXXZ, param::Parameter, obs::MCObservableSet)
    nsites = numsites(model)
    T = convert(Float64, param["T"])
    beta = 1.0 / T

    jk = jackknife(obs)

    for oname in ("Magnetization", "|Magnetization|",
                  "Magnetization^2", "Magnetization^4",
                  "Energy", "Energy^2")
        jk[oname] = jk["Sign * $oname"] / jk["Sign"]
    end

    jk["Binder Ratio"] = jk["Magnetization^4"] / (jk["Magnetization^2"]^2)
    jk["Susceptibility"] = (nsites * beta) * jk["Magnetization^2"]
    jk["Connected Susceptibility"] = (nsites * beta) *
                                     (jk["Magnetization^2"] - jk["|Magnetization|"]^2)
    jk["Specific Heat"] = (nsites * beta * beta) * (jk["Energy^2"] - jk["Energy"]^2)
    return jk
end
