const P = Dict{Symbol, Any}
mutable struct LatticeParameter
    dict :: P
    function LatticeParameter(p::P)
        ld = new(P())
        for key in keys(p)
            val = p[key]
            if isa(val,AbstractArray)
                if isa(first(val), P)
                    setproperty!(ld, key, LatticeParameter.(val))
                else
                    setproperty!(ld, key, val)
                end
            elseif isa(val,P)
                setproperty!(ld, key, LatticeParameter(val))
            else
                setproperty!(ld, key, val)
            end
        end
        return ld
    end
end

import Base: getindex, setindex!, getproperty, setproperty!, propertynames, convert

getproperty(ld::LatticeParameter, name::Symbol) = getfield(ld, :dict)[name]
setproperty!(ld::LatticeParameter, name::Symbol, value) = getfield(ld, :dict)[name] = value
getindex(ld::LatticeParameter, name::Symbol) = ld.name
setindex!(ld::LatticeParameter, name::Symbol, value) = ld.name = value
propertynames(ld::LatticeParameter) = keys(getfield(ld,:dict))

convert(::Type{LatticeParameter}, p::P) = LatticeParameter(p)

