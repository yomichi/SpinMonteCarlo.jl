using JSON
using ChainDiag # https://github.com/yomichi/ChainDiag.jl

function foo(filename, S, T, L, Jz, Jxy, G)
    param = Dict{String,Any}("S" => S,
                             "T" => T,
                             "L" => L,
                             "Jz" => Jz,
                             "Jxy" => Jxy,
                             "Gamma" => G)
    solver = SpinChainSolver(S, L; Jz=Jz, Jxy=Jxy, h=0.0, Guni=G, Gstag=0.0)
    obs = solve(solver, 1.0 / T, 10)
    obs2 = Dict{String,Any}()
    obs2["Energy"] = obs["Energy"]
    obs2["Energy^2"] = obs["Energy^2"]
    obs2["Specific Heat"] = obs["Specific Heat"]
    ## "Order Parameter" and "Susceptibility" are vectors over k = 0:2:L;
    ## the first element (k = 0) is the uniform one measured by the estimator.
    obs2["Magnetization"] = obs["Order Parameter"][1]
    obs2["Susceptibility"] = obs["Susceptibility"][1]
    result = Dict("Parameter" => param, "Result" => obs2)
    open(filename, "w") do io
        return JSON.print(io, result, 2)
    end
end

#                           S,   T, L,   Jz,  Jxy,   G
for (id, p) in enumerate(((0.5, 1.0, 6, 1.0, 1.0, 0.0),
                          (0.5, 1.0, 6, 0.0, 1.0, 0.0),
                          (0.5, 1.0, 6, 1.0, 0.0, 0.0),
                          (0.5, 1.0, 6, -1.0, 0.0, 0.0),
                          (0.5, 1.0, 6, -1.0, -1.0, 0.0),
                          (1.0, 1.0, 6, 1.0, 1.0, 0.0),
                          (0.5, 1.0, 3, 0.0, 1.0, 0.0),
                          (1.5, 1.0, 6, 1.0, 1.0, 0.0),
                          (2.0, 1.0, 4, 1.0, 1.0, 0.0)))
    foo("chain_$id.json", p...)
end
