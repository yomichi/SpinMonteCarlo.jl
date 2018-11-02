using JSON
using ChainDiag # https://github.com/yomichi/ChainDiag.jl

function foo(filename, S,T,L,Jz,Jxy,G)
    param = Dict{String,Any}("S"=>S,
                             "T"=>T,
                             "L"=>L,
                             "Jz"=>Jz,
                             "Jxy"=>Jxy,
                             "Gamma"=>G
                            )
    solver = SpinChainSolver(S,L,Jz=Jz,Jxy=Jxy,h=0.0,Guni=Guni,Gstag=Gstag)
    obs = solve(solver,1.0/T,10)
    obs2 = Dict{String,Any}()
    obs2["Energy"] = obs["Energy"]
    obs2["Energy^2"] = obs["Energy^2"]
    obs2["Specific Heat"] = obs["Specific Heat"]
    obs2["Magnetization"] = obs["Order Parameter"]
    obs2["Susceptibility"] = obs["Susceptibility"]
    result = Dict("Parameter"=>param, "Result"=>obs2)
    open(filename,"w") do io
        JSON.print(io, result, 2)
    end
end

#                           S,   T, L,   Jz,  Jxy,   G
for (id,p) in enumerate(((0.5, 1.0, 6,  1.0,  1.0, 0.0),
                         (0.5, 1.0, 6,  0.0,  1.0, 0.0),
                         (0.5, 1.0, 6,  1.0,  0.0, 0.0),
                         (0.5, 1.0, 6, -1.0,  0.0, 0.0),
                         (0.5, 1.0, 6, -1.0, -1.0, 0.0),
                         (1.0, 1.0, 6,  1.0,  1.0, 0.0),
                         (0.5, 1.0, 3,  0.0,  1.0, 0.0),
                        ))
    foo("chain_$id.json", p...)
end
