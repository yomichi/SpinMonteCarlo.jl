doc"""
Union-find algorithm.
"""
mutable struct UnionFind
    parents :: Vector{Int}
    weights :: Vector{Int}
    ids :: Vector{Int}
    nnodes :: Int
    nclusters :: Int
    UnionFind(n::Int=0) = new(collect(1:n),ones(Int,n),zeros(Int,n),n,n)
end

numnodes(u::UnionFind) = u.nnodes
numclusters(u::UnionFind) = u.nclusters

isroot(u::UnionFind, n::Integer) = u.parents[n] == n

function root_and_weight(u::UnionFind, n::Integer)
    r = n
    @inbounds while !isroot(u,r)
        r = u.parents[r]
    end
    return r,u.weights[r]
end

root(u::UnionFind, n::Integer) = root_and_weight(u,n)[1]

doc"""
    unify!(u, n1, n2)

Connects `n1` and `n2` nodes and returns the root.
"""
function unify!(u::UnionFind, n1::Integer, n2::Integer)
    r1,w1 = root_and_weight(u,n1)
    r2,w2 = root_and_weight(u,n2)
    if r1 != r2
        u.nclusters -= 1
        if w1<w2
            @inbounds u.parents[r1] = r2
            return r2
        elseif w1 == w2
            @inbounds u.parents[r2] = r1
            @inbounds u.weights[r1]+=1
        else
            @inbounds u.parents[r2] = r1
        end
    end
    return r1
end

doc"""
    addnode!(u::UnionFind)

Adds a new node into `u` and returns the number of nodes including the added node.
"""
function addnode!(u::UnionFind) 
    push!(u.parents,length(u.parents)+1)
    push!(u.weights,1)
    push!(u.ids, 0)
    u.nclusters += 1
    u.nnodes += 1
    return u.nnodes
end

doc"""
    clusterize!(u)

Assigns cluster ID to each node and returns the number of clusters.
"""
function clusterize!(u::UnionFind)
    u.nclusters = 0
    @inbounds for i in 1:length(u.parents)
        if isroot(u,i)
            u.nclusters += 1
            u.ids[i] = u.nclusters
        end
    end
    @inbounds for i in 1:length(u.parents)
        u.ids[i] = u.ids[root(u,i)]
    end
    return u.nclusters
end

doc"""
    clusterid(u,i)

Returns the index of the cluster where `i` node belongs.
"""
clusterid(u::UnionFind, i::Integer) = u.ids[i]

