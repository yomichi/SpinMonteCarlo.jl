const Parameter = Dict{String, Any}

doc"""
    convert_parameter(model, param)

generates arguments of updater and estimator from `param` for `model`.

# Example

``` julia-repl
julia> model = Ising(chain_lattice(8));

julia> param = Parameter("T"=>1.0, "J"=>1.0)
Dict{String,Any} with 2 entries:
  "J" => 1.0
  "T" => 1.0

julia> convert_parameter(model, param)
(1.0, [1.0, 1.0])  # `T` and `Js`

julia> param = Parameter("T"=>1.0, "J"=>[1.5,0.5])
Dict{String,Any} with 2 entries:
  "J" => [1.5, 0.5]
  "T" => 1.0

julia> p = convert_parameter(model, param)
(1.0, [1.5, 0.5])

julia> model.spins
8-element Array{Int64,1}:
 -1
 -1
  1
  1
  1
  1
 -1
 -1

julia> local_update!(model, p...);

julia> model.spins
8-element Array{Int64,1}:
 -1
 -1
  1
  1
  1
  1
  1
  1

```

"""
function convert_parameter end

function convert_parameter(model::Union{Ising, Potts, Clock, XY}, param::Parameter)
    T = Float64(param["T"])
    nbt = numbondtypes(model)

    J = get(param, "J", 1.0)
    Js = zeros(nbt)
    Js .= J
    return T, Js
end

function convert_parameter(model::QuantumXXZ, param::Parameter)
    T = Float64(param["T"])

    nbt = numbondtypes(model)
    nst = numsitetypes(model)

    Jz = get(param, "Jz", 1.0)
    Jzs = zeros(nbt)
    Jzs .= Jz

    Jxy = get(param, "Jxy", 1.0)
    Jxys = zeros(nbt)
    Jxys .= Jxy

    G = get(param, "Gamma", 0.0)
    Gs = zeros(nbt)
    Gs .= G

    return T, Jzs, Jxys, Gs
end

