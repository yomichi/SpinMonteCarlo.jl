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
- `clustermag` : Signed magnetization of each cluster before cluster flips.
"""
mutable struct SWInfo
    activated_bonds::Vector{Int}
    clustersize::Vector{Int}
    clusterspin::Vector{Int}
    clustermag::Vector{Int}
end
SWInfo(activated_bonds, clustersize, clusterspin) = SWInfo(activated_bonds, clustersize,
                                                           clusterspin,
                                                           zeros(Int, length(clustersize)))
numclusters(sw::SWInfo) = length(sw.clustersize)

@doc """
    SW_update!(model, param::Parameter)

Updates spin configuration by Swendsen-Wang algorithm

!!! note "Frustrated systems"
    For `Ising`, cluster updates remain exact for any sign of couplings
    (antiferromagnetic and mixed-sign cases included): only satisfied bonds
    (`J σ σ > 0`) are activated, and whole clusters are flipped.
    On unfrustrated (bipartite) lattices an antiferromagnet is gauge-equivalent
    to a ferromagnet, and clusters track the staggered correlations, so the usual
    acceleration is retained.
    On frustrated lattices (loops with an odd number of antiferromagnetic bonds,
    e.g. an antiferromagnet on the triangular lattice), results are still unbiased,
    but clusters decouple from the physical correlations and percolate even at high
    temperature; expect no speedup over `local_update!` (this is a relaxation issue,
    not a correctness issue).
"""
@inline function SW_update!(model::Model, param::Parameter)
    p = convert_parameter(model, param)
    return SW_update!(model, p...)
end

@doc """
    Wolff_update!(model, param::Parameter)

Updates spin configuration by Wolff algorithm

!!! note "Frustrated systems"
    The same caveat as [`SW_update!`](@ref) applies: exact for any sign of
    couplings, but no acceleration should be expected on frustrated lattices.
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
