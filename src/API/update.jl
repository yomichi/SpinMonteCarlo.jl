export local_update!, SW_update!, Wolff_update!, loop_update!

@doc """
    local_update!(model, param)

Updates spin configuration by local spin flip and Metropolice algorithm 
"""
@inline function local_update!(model::Model, param::Parameter)
    p = convert_parameter(model, param)
    return local_update!(model, p...)
end

@doc doc"""
Information of clusters in Swendsen-Wang algorithm.

# Fields
- `activated_bonds` : The number of activated (connected) bonds of each cluster.
- `clustersize` : The number of sites in each cluster.
- `clusterspin` : Spin variable of each cluster (e.g., 1 or -1 for `Ising`).
"""
mutable struct SWInfo
    activated_bonds :: Vector{Int}
    clustersize :: Vector{Int}
    clusterspin :: Vector{Int}
end
numclusters(sw::SWInfo) = length(sw.clustersize)

@doc """
    SW_update!(model, param::Parameter)
    
Updates spin configuration by Swendsen-Wang algorithm
"""
@inline function SW_update!(model::Model, param::Parameter)
    p = convert_parameter(model, param)
    return SW_update!(model, p...)
end

@doc """
    Wolff_update!(model, param::Parameter)

Updates spin configuration by Wolff algorithm
"""
@inline function Wolff_update!(model::Model, param::Parameter)
    p = convert_parameter(model, param)
    return Wolff_update!(model, p...)
end

@doc """
    loop_update!(model, param::Parameter)

Updates spin configuration by loop algorithm 
"""
@inline function loop_update!(model::Model, param::Parameter)
    p = convert_parameter(model, param)
    return loop_update!(model, p...)
end
