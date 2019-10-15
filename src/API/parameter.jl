using MacroTools

export Parameter
export @gen_convert_parameter, convert_parameter

@doc doc"""
    @gen_convert_parameter(model_typename, (keyname, size_fn, default)...)

Generates `convert_parameter(model::model_typename, param::Parameter)`.

# Example

``` julia
@gen_convert_parameter(A, ("A", numbondtypes, 1.0), ("B", 1, 1))
```

generates a function equivalent to the following:

``` julia
doc\"""
    convert_parameter(model::A, param::Parameter)

# Keynames

- "A": a vector with `numbondtypes(model)` elements (default: 1)
- "B": a scalar (default: 1.0)

\"""
function convert_parameter(model::A, param::Parameter)
    ## if `size_fn` is a `Function`,
    ## result is a vector whose size is `size_fn(model)`.
    ## `param["A"]` can take a scalar or a vector.
    a = get(param, "A", 1.0)
    as = zeros(Float64, numbondtypes(model))
    as .= a

    ## otherwise,
    ## result is a scalar.
    b = convert(Int, get(param, "B", 1))

    return as, b
end
```
"""
macro gen_convert_parameter(model_typename, args...)
    body = Expr(:block)
    str_model = sprint(Base.show_unquoted,model_typename)

    document = "    convert_parameter(model::$(str_model), param::Parameter)\n\n"
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
                  esc(:( $(sym) = convert($eltyp, get(param, $name, $default))))
                 )
            push!(syms, sym)
            document *= "- \"$name\": a scalar (default: $default).\n"
        end
    end
    ret = Tuple( s for s in syms )
    push!(body.args, 
          esc(:(return tuple($(syms...))  )) 
         )

    # res_fn = :(function $(esc(:convert_parameter))($(esc(:(model::$model_typename))),
    res_fn = :(function $(esc(:convert_parameter))($(esc(:model))::$(esc(model_typename)),
                                                   $(esc(:param))::$(esc(:Parameter))
                                                  )
                   $body
               end)


    :( @doc $document $res_fn )
end

@doc doc"""
    convert_parameter(model, param)

Generates arguments of updater and estimator.

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
