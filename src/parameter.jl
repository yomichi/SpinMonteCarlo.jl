const Parameter = Dict{String, Any}

using MacroTools

doc"""
    combinedef

similar to `MacroTools.combinedef`, but ignores `rtype`
"""
function combinedef(dict::Dict)
    ## This function is copied and modified from MacroTools.jl
    params = get(dict, :params, [])
    wparams = get(dict, :whereparams, [])
    name = dict[:name]
    name_param = isempty(params) ? name : :($name{$(params...)})
    # We need the `if` to handle parametric inner/outer constructors like
    # SomeType{X}(x::X) where X = SomeType{X}(x, x+2)
    if isempty(wparams)
        :(function $name_param($(dict[:args]...);
                               $(dict[:kwargs]...))
              $(dict[:body])
          end)
    else
        :(function $name_param($(dict[:args]...);
                               $(dict[:kwargs]...)) where {$(wparams...)}
              $(dict[:body])
          end)
    end
end

doc"""
    @gen_convert_parameter(model_typename, (keyname, size_fn, default)...)

generates `convert_parameter(model::model_typename, param::Parameter)`.

# Example

``` julia
@gen_convert_parameter(A, ("A", 1, 1.0), ("B", numbondtypes, 1))
```

generates a function equivalent to the following:

``` julia
doc\"""
    convert_parameter(model::A, param::Parameter)

# Keynames

- "A": a scalar (default: 1.0)
- "B": a vector with `numbondtypes(model)` elements (default: 1)

\"""
function convert_parameter(model::A, param::Parameter)
    ## result is a scalar if `size_fn` isn't `Function`.
    a = get(param, "A", 1.0) :: Float64

    ## result is a vector whose size is `size_fn(model)`.
    ## `param["B"]` can take a scalar or a vector.
    b = get(param, "B", 1)
    bs = zeros(Int, numbondtypes(model))
    bs .= b

    return a, bs
end
```
"""
macro gen_convert_parameter(model_typename, args...)
    fn = Dict{Symbol, Any}()
    fn[:name] = esc(:convert_parameter)
    fn[:args] = Any[esc(:(model::$model_typename)), esc(:(param::Parameter))]
    fn[:kwargs] = []
    fn[:whereparams] = ()

    body = Expr(:block)

    document = "    convert_parameter(model::$(model_typename), param::Parameter)\n\n"
    document *= "# Keynames:\n"

    syms = Symbol[]
    for arg in args
        @capture(arg, (name_String, sz_, default_))
        eltyp = typeof(default)
        sym = gensym()
        if isa(eval(sz), Function)
            push!(body.args,
                  esc(:( $(sym) = get(param, $name, $default) ))
                 )
            sym2 = gensym()
            push!(syms, sym2)
            push!(body.args,
                  esc(:( $(sym2) = zeros($eltyp, $sz(model))))
                 )
            push!(body.args,
                  esc(:( $(sym2) .= $(sym) ))
                 )
            document *= "- \"$name\": a vector with `$(sz)(model)` elements (default: $default).\n"
        else
            push!(body.args,
                  esc(:( $(sym) = get(param, $name, $default) :: $eltyp ))
                 )
            push!(syms, sym)
            document *= "- \"$name\": a scalar (default: $default).\n"
        end
    end
    ret = Tuple( s for s in syms )
    push!(body.args, 
          # esc(:( return $(syms...) ))
          esc(:(return tuple($(syms...))  )) 
         )
    fn[:body] = body

    res_fn = combinedef(fn)
    :( @doc $document $res_fn )
end

doc"""
    convert_parameter(model, param)

generates arguments of updater and estimator.

# Example
``` julia-repl
julia> model = Ising(chain_lattice(4));

julia> p = convert_parameter(model, Parameter("J"=>1.0))
(1.0, [1.0, 1.0]) # T and Js

julia> p = convert_parameter(model, Parameter("J"=>[1.5, 0.5]))
(1.0, [1.5, 0.5]) # J can take a vector whose size is `numbondtypes(model)`

julia> model.spins
4-element Array{Int64,1}:
  1
  1
  1
 -1

julia> local_update!(model, p...);

julia> model.spins
4-element Array{Int64,1}:
  1
  1
 -1
 -1
```
"""
function convert_parameter end

@gen_convert_parameter(Union{Ising, Potts, Clock, XY}, ("T", 1, 1.0),
                                                       ("J", numbondtypes, 1.0),
                                                      )
@gen_convert_parameter(QuantumXXZ, ("T", 1, 1.0),
                                   ("Jz", numbondtypes, 1.0),
                                   ("Jxy", numbondtypes, 1.0),
                                   ("Gamma", numsitetypes, 0.0),
                                  )
