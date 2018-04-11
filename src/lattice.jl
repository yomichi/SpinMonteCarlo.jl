"""
    Lattice

Representing a lattice structure.

# Arguments

- `dim :: Int`: Dimension of the lattice
- `size :: Vector{Int}`: Size
- `nsitetypes :: Int`: The number of site types
- `nbondtypes :: Int`: The number of bond types
- `sites :: Vector{Vector{Int}}`: Site indices belonging each site type
- `bonds :: Vector{Vector{Int}}`: Bond indices belonging each bond type
- `nsites :: Int`: The number of sites
- `nbonds :: Int`: The number of bonds
- `sitetypes :: Vector{Int}`: Site type of each site
- `bondtypes :: Vector{Int}`: Bond type of each bond
- `transvector :: Matrix{Float64}`: Translational vectors as a `dim` x `dim` matrix
- `site_coords :: Matrix{Float64}`: Coordinate of each site  as a `dim` x `nsites` matrix
- `bond_dirs :: Matrix{Float64}`: Direction of each bond as a `dim` x `nbonds` matrix
- `neighborsites :: Vector{Vector{Int}}`: Indices of neighbor sites of each site
- `neighborbonds :: Vector{Vector{Int}}`: Indices of adjacent bonds of each site
- `source :: Vector{Int}`: Index of a site connected to each bond
- `target :: Vector{Int}`: Index of another site connected to each bond
- `site_L2 :: Vector{Int}`: Index of the site that is `size/2` from each site
- `site_L4 :: Vector{Int}`: Index of the site that is `size/4` from each site

"""
mutable struct Lattice
    "dim"
    dim :: Int
    size :: Vector{Int}
    nsitetypes :: Int
    nbondtypes :: Int
    sites :: Vector{Vector{Int}}
    bonds :: Vector{Vector{Int}}
    nsites :: Int
    nbonds :: Int
    sitetypes :: Vector{Int}
    bondtypes :: Vector{Int}
    transvector :: Matrix{Float64}
    site_coords :: Matrix{Float64}
    bond_dirs :: Matrix{Float64}
    neighborsites :: Vector{Vector{Int}}
    neighborbonds :: Vector{Vector{Int}}
    source :: Vector{Int}
    target :: Vector{Int}
    site_L2 :: Vector{Int}
    site_L4 :: Vector{Int}
end

import Base.size
import Distributions.dim

"""
    dim(lat::Lattice)
    dim(model::Model)

return the dimension of lattice.
"""
dim(lat::Lattice) = lat.dim
dim(model::Model) = dim(model.lat)

"""
    size(lat::Lattice, [dim::Integer])
    size(model::Model, [dim::Integer])

return the size of lattice.
"""
size(lat::Lattice) = lat.size
size(lat::Lattice, dim::Integer) = lat.size[dim]
size(model::Model) = size(model.lat)
size(model::Model, dim::Integer) = size(model.lat,dim)

"""
    sites(lat::Lattice, sitetype::Integer)
    sites(model::Model, sitetype::Integer)

return sites with `sitetype`
"""
sites(lat::Lattice, sitetype::Integer) = lat.sites[sitetype]
sites(model::Model, sitetype::Integer) = sites(model.lat, sitetype)

"""
    numsites(lat::Lattice)
    numsites(model::Model)

return the number of all sites.
"""
numsites(lat::Lattice) = lat.nsites
numsites(model::Model) = numsites(model.lat)

"""
    numsites(lat::Lattice, sitetype::Integer)
    numsites(model::Model, sitetype::Integer)

return the number of `sitetype` sites.
"""
numsites(lat::Lattice, sitetype::Integer) = length(lat.sites[sitetype])
numsites(model::Model, sitetype::Integer) = numsites(model.lat, sitetype)

"""
    bonds(lat::Lattice, bondtype::Integer)
    bonds(model::Model, bondtype::Integer)

return bonds with `bondtype`
"""
bonds(lat::Lattice, bondtype::Integer) = lat.bonds[bondtype]
bonds(model::Model, bondtype::Integer) = bonds(model.lat, bondtype)

"""
    numbonds(lat::Lattice)
    numbonds(model::Model)

return the number of all bonds.
"""
numbonds(lat::Lattice) = lat.nbonds
numbonds(model::Model) = numbonds(model.lat)

"""
    numbonds(lat::Lattice, bondtype::Integer)
    numbonds(model::Model, bondtype::Integer)

return the number of `bondtype` bonds.
"""
numbonds(lat::Lattice, bondtype::Integer) = length(lat.bonds[bondtype])
numbonds(model::Model, bondtype::Integer) = numbonds(model.lat, bondtype)

"""
    numsitetypes(lat::Lattice)
    numsitetypes(model::Model)

return the number of sitetypes.
"""
numsitetypes(lat::Lattice) = lat.nsitetypes
numsitetypes(model::Model) = numsitetypes(model.lat)

"""
    numbondtypes(lat::Lattice)
    numbondtypes(model::Model)

return the number of bondtypes.
"""
numbondtypes(lat::Lattice) = lat.nbondtypes
numbondtypes(model::Model) = numbondtypes(model.lat)

"""
    neighborsites(lat::Lattice, site::Integer)
    neighborsites(model::Model, site::Integer)

return the neighbor sites of `site`.
"""
neighborsites(lat::Lattice, site::Integer) = lat.neighborsites[site]
neighborsites(model::Model, site::Integer) = neighborsites(model.lat, site)

"""
    neighborbonds(lat::Lattice, site::Integer)
    neighborbonds(model::Model, site::Integer)

return the neighbor bonds of `site`.
"""
neighborbonds(lat::Lattice, site::Integer) = lat.neighborbonds[site]
neighborbonds(model::Model, site::Integer) = neighborbonds(model.lat, site)

"""
    neighbors(lat::Lattice, site::Integer)
    neighbors(model::Model, site::Integer)

return the neighbor sites and bonds of `site`.
"""
neighbors(lat::Lattice, site::Integer) = zip(neighborsites(lat,site), neighborbonds(lat,site))
neighbors(model::Model, site::Integer) = neighbors(model.lat, site)

"""
    source(lat::Lattice, bond::Integer)
    source(model::Model, bond::Integer)

return the source site of `bond`.
"""
source(lat::Lattice, bond::Integer) = lat.source[bond]
source(model::Model, bond::Integer) = source(model.lat, bond)

"""
    target(lat::Lattice, bond::Integer)
    target(model::Model, bond::Integer)

return the target site of `bond`.
"""
target(lat::Lattice, bond::Integer) = lat.target[bond]
target(model::Model, bond::Integer) = target(model.lat, bond)

"""
    sitecoordinate(lat::Lattice, site::Integer)
    sitecoordinate(model::Model, site::Integer)

return the coordinate of the `site` in the Cartesian system
"""
sitecoordinate(lat::Lattice, site::Integer) = lat.transvector * lat.site_coords[:,site]
sitecoordinate(model::Model, site::Integer) = sitecoordinate(model.lat, site)

"""
    lattice_sitecoordinate(lat::Lattice, site::Integer)
    lattice_sitecoordinate(model::Model, site::Integer)

return the coordinate of the `site` in the lattice system
"""
lattice_sitecoordinate(lat::Lattice, site::Integer) = lat.site_coords[:, site]
lattice_sitecoordinate(model::Model, site::Integer) = lattice_sitecoordinate(model.lat, site)

"""
    sitetype(lat::Lattice, site::Integer)
    sitetype(model::Model, site::Integer)

return the type of `site`
"""
sitetype(lat::Lattice, site::Integer) = lat.sitetypes[site]
sitetype(model::Model, site::Integer) = sitetype(model.lat, site)

"""
    bondtype(lat::Lattice, bond::Integer)
    bondtype(model::Model, bond::Integer)

return the type of `bond`
"""
bondtype(lat::Lattice, bond::Integer) = lat.bondtypes[bond]
bondtype(model::Model, bond::Integer) = bondtype(model.lat, bond)

"""
    bonddirection(lat::Lattice, bond::Integer)
    bonddirection(model::Model, bond::Integer)

return the direction of the `bond` as vector in the Cartesian system
"""
bonddirection(lat::Lattice, bond::Integer) = lat.transvector * lat.bond_dirs[:, bond]
bonddirection(model::Model, bond::Integer) = bonddirection(model.lat, bond)

"""
    lattice_bonddirection(lat::Lattice, bond::Integer)
    lattice_bonddirection(model::Model, bond::Integer)

return the direction of the `bond` as vector in the lattice system
"""
lattice_bonddirection(lat::Lattice, bond::Integer) = lat.bond_dirs[:, bond]
lattice_bonddirection(model::Model, bond::Integer) = lattice_bonddirection(model.lat, bond)

"""
    siteL2(lat::Lattice, site::Integer)
    siteL2(model::Model, site::Integer)

return (x+L/2, y+W/2, ...) site
"""
siteL2(lat::Lattice, site::Integer) = lat.site_L2[site]
siteL2(model::Model, site::Integer) = siteL2(model.lat, site)

"""
    siteL4(lat::Lattice, site::Integer)
    siteL4(model::Model, site::Integer)

return (x+L/4, y+W/4, ...) site
"""
siteL4(lat::Lattice, site::Integer) = lat.site_L4[site]
siteL4(model::Model, site::Integer) = siteL4(model.lat, site)

"""
    dimer_lattice()
    dimer_lattice(params::Dict)
    
generate dimer lattice (indeed, this is not a "lattice")
"""
function dimer_lattice()
    dim = 1
    sz = [2]
    nsitetypes = 2
    nbondtypes = 1
    sites = Vector{Int}[zeros(Int,0) for i in 1:nsitetypes]
    bonds = Vector{Int}[zeros(Int,0) for i in 1:nbondtypes]
    nsites = 2
    nbonds = 1
    sitetypes = [1,2]
    bondtypes = [1]
    transvector = eye(dim)
    coords = zeros(dim, 2)
    bond_dirs = zeros(dim, 1)
    neighborsites = [[2],[1]]
    neighborbonds = [[1],[1]]
    source = zeros(Int,1)
    target = zeros(Int,1)
    site_L2 = zeros(Int,2)
    site_L4 = zeros(Int,2)

    sites[1] = [1]
    sites[2] = [2]
    bonds[1] = [1]
    coords[1,1] = 0.0
    coords[1,2] = 1.0
    bond_dirs[1,1] = 1.0
    source[1] = 1
    target[1] = 2
    site_L2[1] = 2
    site_L2[2] = 1
    site_L4[1] = 2
    site_L4[2] = 1

    Lattice(dim, sz, nsitetypes, nbondtypes,
            sites, bonds, nsites, nbonds, sitetypes, bondtypes,
            transvector, coords, bond_dirs,
            neighborsites, neighborbonds,
            source,target, site_L2, site_L4)
end
dimer_lattice(params::Dict) = dimer_lattice()

"""
    chain_lattice(L::Integer)
    chain_lattice(params::Dict)
    
generate chain lattice with length `L` or params["L"].
"""
function chain_lattice(L::Integer)
    dim = 1
    sz = [L]
    nsitetypes = 2
    nbondtypes = 2
    nsites = L
    nbonds = L
    sitetypes = zeros(Int,L)
    bondtypes = zeros(Int,L)
    transvector = eye(dim)
    coords = zeros(dim, L)
    bond_dirs = zeros(dim, L)
    neighborsites = Vector{Int}[]
    neighborbonds = Vector{Int}[]
    source = zeros(Int,L)
    target = zeros(Int,L)
    site_L2 = zeros(Int,L)
    site_L4 = zeros(Int,L)

    ns = zeros(Int,2)
    nb = zeros(Int,2)
    for s in 1:L
        coords[1,s] = (s-1)
        bond_dirs[1,s] = 1.0
        ns[1] = mod1(s+1,L)
        ns[2] = mod1(s-1,L)
        nb[1] = s
        nb[2] = mod1(s-1,L)
        push!(neighborsites, copy(ns))
        push!(neighborbonds, copy(nb))
        source[s] = s
        target[s] = neighborsites[s][1]
        sitetypes[s] = ifelse(isodd(s),1,2)
        bondtypes[s] = ifelse(isodd(s),1,2)
        site_L2[s] = mod1(s+div(L,2),L)
        site_L4[s] = mod1(s+div(L,4),L)
    end
    sites = Vector{Int}[zeros(Int,0) for i in 1:nsitetypes]
    bonds = Vector{Int}[zeros(Int,0) for i in 1:nbondtypes]
    for s in 1:nsites
        push!(sites[sitetypes[s]], s)
    end
    for b in 1:nbonds
        push!(bonds[bondtypes[b]], b)
    end
    Lattice(dim, sz, nsitetypes, nbondtypes,
            sites, bonds, nsites, nbonds, sitetypes, bondtypes,
            transvector, coords, bond_dirs, neighborsites, neighborbonds,
            source,target, site_L2, site_L4)
end
chain_lattice(params::Dict) = chain_lattice(params["L"])

"""
    square_lattice(L::Integer, W::Integer=L)
    square_lattice(params::Dict)
    
generate square lattice with size `L` \\times `W`.
"""
square_lattice(L::Integer) = square_lattice(L,L)
function square_lattice(L::Integer, W::Integer)
    s2xy(s::Integer) = mod(s-1,L),div(s-1,L)
    xy2s(x::Integer, y::Integer) = mod(y,W)*L + mod(x,L) + 1

    dim = 2
    sz = [L,W]
    nsitetypes = 2
    nbondtypes = 2
    nsites = L*W
    nbonds = 2*L*W
    sitetypes = zeros(Int, nsites)
    bondtypes = zeros(Int, nbonds)
    transvector = eye(dim)
    coords = zeros(dim, nsites)
    bond_dirs = zeros(dim, nbonds)
    neighborsites = Vector{Int}[]
    neighborbonds = Vector{Int}[]
    source = zeros(Int,nbonds)
    target = zeros(Int,nbonds)
    sitetypes = zeros(Int,nsites)
    bondtypes = zeros(Int,nbonds)
    site_L2 = zeros(Int,nsites)
    site_L4 = zeros(Int,nsites)
    L2 = div(L,2)
    L4 = div(L,4)
    W2 = div(W,2)
    W4 = div(W,4)
    ns = zeros(Int,4)
    nb = zeros(Int,4)
    @inbounds for s in 1:nsites
        x,y = s2xy(s)
        ns[1] = xy2s(x+1,y)
        ns[2] = xy2s(x,y+1)
        ns[3] = xy2s(x-1,y)
        ns[4] = xy2s(x,y-1)
        push!(neighborsites, copy(ns))
        nb[1] = 2s-1
        nb[2] = 2s
        nb[3] = 2ns[3]-1
        nb[4] = 2ns[4]
        push!(neighborbonds, copy(nb))
        source[2s-1:2s] .= s
        target[2s-1:2s] .= neighborsites[s][1:2]
        coords[1,s] = x
        coords[2,s] = y
        sitetypes[s] = ifelse(iseven(x+y),1,2)
        bond_dirs[:, 2s-1] = [1.0, 0.0]
        bond_dirs[:, 2s] = [0.0, 1.0]
        bondtypes[2s-1] = 1 
        bondtypes[2s] = 2 
        site_L2[s] = xy2s(x+L2, y+W2)
        site_L4[s] = xy2s(x+L4, y+W4)
    end
    sites = Vector{Int}[zeros(Int,0) for i in 1:nsitetypes]
    bonds = Vector{Int}[zeros(Int,0) for i in 1:nbondtypes]
    for s in 1:nsites
        push!(sites[sitetypes[s]], s)
    end
    for b in 1:nbonds
        push!(bonds[bondtypes[b]], b)
    end
    Lattice(dim, sz, nsitetypes, nbondtypes,
            sites, bonds, nsites, nbonds, sitetypes, bondtypes,
            transvector, coords, bond_dirs, neighborsites, neighborbonds,
            source,target, site_L2, site_L4)
end
function square_lattice(params::Dict)
    L = params["L"]
    W = get(params, "W", L)
    return square_lattice(L, W)
end

"""
    triangular_lattice(L::Integer, W::Integer=L)
    triangular_lattice(params::Dict)
    
generate triangular lattice with size `L` \\times `W`.
"""
triangular_lattice(L::Integer) = triangular_lattice(L,L)
function triangular_lattice(L::Integer, W::Integer)
    s2xy(s::Integer) = mod(s-1,L),div(s-1,L)
    xy2s(x::Integer, y::Integer) = mod(y,W)*L + mod(x,L) + 1
    ex = [1.0, 0.0]
    ey = [0.5, 0.5sqrt(3)]

    dim = 2
    sz = [L,W]
    nsitetypes = 3
    nbondtypes = 3
    nsites = L*W
    nbonds = 3nsites
    sitetypes = zeros(Int,nsites)
    bondtypes = zeros(Int,nbonds)
    transvector = [1.0 0.5; 0.0 0.5sqrt(3)]
    coords = zeros(dim, nsites)
    bond_dirs = zeros(dim, nbonds)
    neighborsites = Vector{Int}[]
    neighborbonds = Vector{Int}[]
    source = zeros(Int,nbonds)
    target = zeros(Int,nbonds)
    site_L2 = zeros(Int, nsites)
    site_L4 = zeros(Int, nsites)
    L2 = div(L,2)
    L4 = div(L,4)
    W2 = div(W,2)
    W4 = div(W,4)
    ns = zeros(Int,6)
    nb = zeros(Int,6)
    @inbounds for s in 1:nsites
        x,y = s2xy(s)
        ns[1] = xy2s(x+1,y)
        ns[2] = xy2s(x,y+1)
        ns[3] = xy2s(x-1,y+1)
        ns[4] = xy2s(x-1,y)
        ns[5] = xy2s(x,y-1)
        ns[6] = xy2s(x+1,y-1)
        push!(neighborsites, copy(ns))
        nb[1] = 3s-2
        nb[2] = 3s-1
        nb[3] = 3s
        nb[4] = 3ns[4]-2
        nb[5] = 3ns[5]-1
        nb[6] = 3ns[6]
        push!(neighborbonds, copy(nb))
        source[3s-2:3s] .= s
        target[3s-2:3s] .= neighborsites[s][1:3]
        coords[1,s] = x
        coords[2,s] = y
        if sitetypes[s]==0
            sitetypes[s] = 1
        else
            if sitetypes[neighbors[1,s]]==0
                sitetypes[neighbors[1,s]]=mod1(sitetypes[s]+1,3)
            end
            if sitetypes[neighbors[2,s]]==0
                sitetypes[neighbors[2,s]]=mod1(sitetypes[s]+2,3)
            end
        end
        bond_dirs[:, 3s-2] .= [1.0,0.0]
        bond_dirs[:, 3s-1] .= [0.0,1.0]
        bond_dirs[:, 3s] .= [-1.0,1.0]
        bondtypes[3s-2] = 1
        bondtypes[3s-1] = 2
        bondtypes[3s] = 3
        site_L2[s] = xy2s(x+L2, y+W2)
        site_L4[s] = xy2s(x+L4, y+W4)
    end
    sites = Vector{Int}[zeros(Int,0) for i in 1:nsitetypes]
    bonds = Vector{Int}[zeros(Int,0) for i in 1:nbondtypes]
    for s in 1:nsites
        push!(sites[sitetypes[s]], s)
    end
    for b in 1:nbonds
        push!(bonds[bondtypes[b]], b)
    end
    Lattice(dim, sz, nsitetypes, nbondtypes,
            sites, bonds, nsites, nbonds, sitetypes, bondtypes,
            transvector, coords, bond_dirs, neighborsites, neighborbonds,
            source,target, site_L2, site_L4)
end
function triangular_lattice(params::Dict)
    L = params["L"]
    W = get(params, "W", L)
    return triangular_lattice(L, W)
end

"""
    cubic_lattice(L::Integer, W::Integer=L, H::Integer=W)
    
generate cubic lattice with size `L` \\times `W` \\times `H`.
"""
cubic_lattice(L::Integer) = cubic_lattice(L,L,L)
cubic_lattice(L::Integer, W::Integer) = cubic_lattice(L,W,W)
function cubic_lattice(L::Integer, W::Integer, H::Integer)
    function s2xyz(s::Integer)
        z = div(s-1, L*W)
        xy = mod(s-1, L*W)
        x = mod(xy, L)
        y = div(xy, L)
        return x, y, z
    end
    xyz2s(x::Integer, y::Integer, z::Integer) = mod(z,H)*W*L + mod(y,W)*L + mod(x,L) + 1

    dim = 3
    sz = [L,W,H]
    nsitetypes = 2
    nbondtypes = 3
    nsites = L*W*H
    nbonds = 3nsites
    sitetypes = zeros(Int,nsites)
    bondtypes = zeros(Int,nbonds)
    transvector = eye(dim)
    coords = zeros(dim, nsites)
    bond_dirs = zeros(dim, nbonds)
    neighborsites = Vector{Int}[]
    neighborbonds = Vector{Int}[]
    source = zeros(Int,nbonds)
    target = zeros(Int,nbonds)
    site_L2 = zeros(Int,nsites)
    site_L4 = zeros(Int,nsites)
    L2 = div(L,2)
    L4 = div(L,4)
    W2 = div(W,2)
    W4 = div(W,4)
    H2 = div(H,2)
    H4 = div(H,4)
    ns = zeros(Int,6)
    nb = zeros(Int,6)
    @inbounds for s in 1:nsites
        x,y,z = s2xyz(s)
        ns[1] = xyz2s(x+1,y,z)
        ns[2] = xyz2s(x,y+1,z)
        ns[3] = xyz2s(x,y,z+1)
        ns[4] = xyz2s(x-1,y,z)
        ns[5] = xyz2s(x,y-1,z)
        ns[6] = xyz2s(x,y,z-1)
        push!(neighborsites, copy(ns))
        nb[1] = 3s-2
        nb[2] = 3s-1
        nb[3] = 3s
        nb[4] = 3ns[4]-2
        nb[5] = 3ns[5]-1
        nb[6] = 3ns[6]
        push!(neighborbonds, copy(nb))
        source[3s-2:3s] .= s
        target[3s-2:3s] .= neighborsites[s][1:3]
        coords[1,s] = x
        coords[2,s] = y
        coords[3,s] = z
        sitetypes[s] = ifelse(iseven(x+y+z),1,2)
        bondtypes[3s-2] = 1
        bondtypes[3s-1] = 2
        bondtypes[3s] = 3
        bond_dirs[:, 3s-2] = [1.0, 0.0, 0.0]
        bond_dirs[:, 3s-1] = [0.0, 1.0, 0.0]
        bond_dirs[:, 3s-0] = [0.0, 0.0, 1.0]
        site_L2[s] = xyz2s(x+L2,y+W2,z+H2)
        site_L4[s] = xyz2s(x+L4,y+W4,z+H4)
    end
    sites = Vector{Int}[zeros(Int,0) for i in 1:nsitetypes]
    bonds = Vector{Int}[zeros(Int,0) for i in 1:nbondtypes]
    for s in 1:nsites
        push!(sites[sitetypes[s]], s)
    end
    for b in 1:nbonds
        push!(bonds[bondtypes[b]], b)
    end
    Lattice(dim, sz, nsitetypes, nbondtypes,
            sites, bonds, nsites, nbonds, sitetypes, bondtypes,
            transvector, coords, bond_dirs,
            neighborsites, neighborbonds,
            source,target, site_L2, site_L4)
end
function cubic_lattice(params::Dict)
    L = params["L"]
    W = get(params, "W", L)
    H = get(params, "H", W)
    return cubic_lattice(L, W, H)
end

"""
    fully_connected_lattice(N::Integer)
    
generate `N` site fully connected lattice
"""
function fully_connected_lattice(N::Integer)
    dim = 1
    nbonds = div(N*(N-1),2)
    nsitetypes = 1
    nbondtypes = 1
    sites = [collect(1:N)]
    bonds = [collect(1:nbonds)]
    sitetypes = ones(N)
    bondtypes = ones(nbonds)
    neighborsites = [vcat(1:s-1, s+1:N) for s in 1:N]
    neighborbonds = Vector{Int}[zeros(Int,0) for s in 1:N]
    source = zeros(Int,nbonds)
    target = zeros(Int,nbonds)
    ib = 0
    for i in 1:N
        for j in (i+1):N
            ib += 1
            source[ib] = i
            target[ib] = j
            push!(neighborbonds[i], ib)
            push!(neighborbonds[j], ib)
        end
    end

    lat = Lattice(dim, [N], nsitetypes, nbondtypes, sites, bonds, N, nbonds, sitetypes, bondtypes,
            zeros(1,1), zeros(1,N), zeros(1,nbonds), neighborsites, neighborbonds, source, target, zeros(Int,N), zeros(Int,N))
    return lat
end
fully_connected_lattice(params::Dict) = fully_connected_lattice(params["N"])
