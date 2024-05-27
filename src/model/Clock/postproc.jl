@doc raw"""
    postproc(model::Clock, param::Parameter, obs::MCObservableSet)

# Observables to be calculated
In the following, $m$ is total magnetization per site and $\epsilon$ is total energy per site.

- `"Binder Ratio x"`
    - $\frac{\left \langle m_x^4 \right \rangle}{\left \langle m_x^2 \right\rangle^2}$
- `"Binder Ratio y"`
    - $\frac{\left \langle m_y^4 \right \rangle}{\left \langle m_y^2 \right\rangle^2}$
- `"Binder Ratio"`
    - $\frac{\left \langle |m|^4 \right \rangle}{\left \langle |m|^2 \right\rangle^2}$
- `"Susceptibility x"`
    - $\frac{N}{T}\left(\left\langle m_x^2\right\rangle\right)$
- `"Susceptibility y"`
    - $\frac{N}{T}\left(\left\langle m_y^2\right\rangle\right)$
- `"Susceptibility y"`
    - $\frac{N}{T}\left(\left\langle |m|^2\right\rangle\right)$
- `"Connected Susceptibility x"`
    - $\frac{N}{T}\left(\left\langle m_x^2\right\rangle - \left\langle |m_x| \right\rangle^2\right)$
- `"Connected Susceptibility y"`
    - $\frac{N}{T}\left(\left\langle m_y^2\right\rangle - \left\langle |m_y| \right\rangle^2\right)$
- `"Connected Susceptibility"`
    - $\frac{N}{T}\left(\left\langle |m|^2\right\rangle - \left\langle |m| \right\rangle^2\right)$
- `"Specific Heat"`
    - $\frac{N}{T^2}\left(\left\langle \epsilon^2\right\rangle - \left\langle \epsilon \right\rangle^2\right)$
"""
function postproc(model::Clock, param::Parameter, obs::MCObservableSet)
    nsites = numsites(model)
    T = convert(Float64, param["T"])
    beta = 1.0 / T

    jk = jackknife(obs)
    jk["Binder Ratio x"] = jk["Magnetization x^4"] / (jk["Magnetization x^2"]^2)
    jk["Binder Ratio y"] = jk["Magnetization y^4"] / (jk["Magnetization y^2"]^2)
    jk["Binder Ratio"] = jk["|Magnetization|^4"] / (jk["|Magnetization|^2"]^2)
    jk["Susceptibility x"] = (nsites * beta) * jk["Magnetization x^2"]
    jk["Susceptibility y"] = (nsites * beta) * jk["Magnetization y^2"]
    jk["Susceptibility"] = (nsites * beta) * jk["|Magnetization|^2"]
    jk["Connected Susceptibility x"] = (nsites * beta) *
                                       (jk["Magnetization x^2"] - jk["|Magnetization x|"]^2)
    jk["Connected Susceptibility y"] = (nsites * beta) *
                                       (jk["Magnetization y^2"] - jk["|Magnetization y|"]^2)
    jk["Connected Susceptibility"] = (nsites * beta) *
                                     (jk["|Magnetization|^2"] - jk["|Magnetization|"]^2)
    jk["Specific Heat"] = (nsites * beta * beta) * (jk["Energy^2"] - jk["Energy"]^2)
    return jk
end
