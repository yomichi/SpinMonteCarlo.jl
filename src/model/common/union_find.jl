export UnionFind, addnode!, unify!, clusterize!, clusterid

@doc doc"""
Union-find algorithm.
"""
mutable struct UnionFind
    parents::Vector{Int}
    weights::Vector{Int}
    ids::Vector{Int}
    nnodes::Int
    nclusters::Int
    UnionFind(n::Int=0) = new(collect(1:n), ones(Int, n), zeros(Int, n), n, n)
end

numnodes(u::UnionFind) = u.nnodes
numclusters(u::UnionFind) = u.nclusters

isroot(u::UnionFind, n::Integer) = u.parents[n] == n

@doc doc"""
    root_path_halving!(u::UnionFind, n::Integer)

Returns the root node of the cluster where `n` belongs.
This may changes graph connection by "Path halving" method.
"""
function root_path_halving!(u::UnionFind, n::Integer)
    @inbounds while !isroot(u, n)
        n = u.parents[n] = u.parents[u.parents[n]]
    end
    return n
end

@doc doc"""
    root_path_splitting!(u::UnionFind, n::Integer)

Returns the root node of the cluster where `n` belongs.
This may changes graph connection by "Path splitting" method.
"""
function root_path_splitting!(u::UnionFind, n::Integer)
    @inbounds while !isroot(u, n)
        p = u.parents[n]
        u.parents[n] = u.parents[p]
        n = p
    end
    return n
end

@doc doc"""
    root!(u::UnionFind, n::Integer)

Returns the root node of the cluster where `n` belongs.
This may changes graph connection by "Path halving" method.
"""
const root! = root_path_halving!

@doc doc"""
    unify!(u, n1, n2)

Connects `n1` and `n2` nodes using union by weight and returns the root.
"""
function unify!(u::UnionFind, n1::Integer, n2::Integer)
    r1 = root!(u, n1)
    r2 = root!(u, n2)
    @inbounds if r1 != r2
        u.nclusters -= 1
        if u.weights[r1] < u.weights[r2]
            r1, r2 = r2, r1
        end
        u.parents[r2] = r1
        u.weights[r1] += u.weights[r2]
    end
    return r1
end

@doc doc"""
    addnode!(u::UnionFind)

Adds a new node into `u` and returns the number of nodes including the added node.
"""
function addnode!(u::UnionFind)
    push!(u.parents, length(u.parents) + 1)
    push!(u.weights, 1)
    push!(u.ids, 0)
    u.nclusters += 1
    u.nnodes += 1
    return u.nnodes
end

@doc doc"""
    clusterize!(u::UnionFind)

Assigns cluster ID to each node and returns the number of clusters.
"""
function clusterize!(u::UnionFind)
    u.nclusters = 0
    @inbounds for i in 1:length(u.parents)
        if isroot(u, i)
            u.nclusters += 1
            u.ids[i] = u.nclusters
        end
    end
    @inbounds for i in 1:length(u.parents)
        u.ids[i] = u.ids[root!(u, i)]
    end
    return u.nclusters
end

@doc doc"""
    clusterid(u::UnionFind, i::Integer)

Returns the index of the cluster where `i` node belongs.
"""
clusterid(u::UnionFind, i::Integer) = u.ids[i]
