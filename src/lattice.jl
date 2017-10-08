type Lattice
    dim :: Int
    size :: Vector{Int}
    nsitetypes :: Int
    nbondtypes :: Int
    nsites :: Vector{Int}
    nbonds :: Vector{Int}
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

return the dimension of lattice.
"""
dim(lat::Lattice) = lat.dim

"""
    size(lat::Lattice, [dim::Integer])

return the size of lattice.
"""
size(lat::Lattice) = lat.size
size(lat::Lattice, dim::Integer) = lat.size[dim]

"""
    numsites(lat::Lattice)

return the number of all sites.
"""
numsites(lat::Lattice) = sum(lat.nsites)

"""
    numsites(lat::Lattice, sitetype::Integer)

return the number of `sitetype` sites.
"""
numsites(lat::Lattice, sitetype::Integer) = lat.nsites[sitetype]

"""
    numbonds(lat::Lattice)

return the number of all bonds.
"""
numbonds(lat::Lattice) = sum(lat.nbonds)

"""
    numbonds(lat::Lattice, bondtype::Integer)

return the number of `bondtype` bonds.
"""
numbonds(lat::Lattice, bondtype::Integer) = lat.nbonds[bondtype]

"""
    numsitetypes(lat::Lattice)

return the number of sitetypes.
"""
numsitetypes(lat::Lattice) = lat.nsitetypes

"""
    numbondtypes(lat::Lattice)

return the number of bondtypes.
"""
numbondtypes(lat::Lattice) = lat.nbondtypes

"""
    neighborsites(lat::Lattice, site::Integer)

return the neighbor sites of `site`.
"""
neighborsites(lat::Lattice, site::Integer) = lat.neighborsites[site]

"""
    neighborbonds(lat::Lattice, site::Integer)

return the neighbor bonds of `site`.
"""
neighborbonds(lat::Lattice, site::Integer) = lat.neighborbonds[site]

"""
    neighbors(lat::Lattice, site::Integer)

return the neighbor sites and bonds of `site`.
"""
neighbors(lat::Lattice, site::Integer) = zip(neighborsites(lat,site), neighborbonds(lat,site))

"""
    source(lat::Lattice, bond::Integer)

return the source site of `bond`.
"""
source(lat::Lattice, bond::Integer) = lat.source[bond]

"""
    target(lat::Lattice, bond::Integer)

return the target site of `bond`.
"""
target(lat::Lattice, bond::Integer) = lat.target[bond]

"""
    sitecoordinate(lat::Lattice, site::Integer)

return the coordinate of the `site` in the Cartesian system
"""
sitecoordinate(lat::Lattice, site::Integer) = lat.transvector * lat.site_coords[:,site]

"""
    lattice_sitecoordinate(lat::Lattice, site::Integer)

return the coordinate of the `site` in the lattice system
"""
lattice_sitecoordinate(lat::Lattice, site::Integer) = lat.site_coords[:, site]

"""
    sitetype(lat::Lattice, site::Integer)
return the type of `site`
"""
sitetype(lat::Lattice, site::Integer) = lat.sitetypes[site]

"""
    bondtype(lat::Lattice, bond::Integer)
return the type of `bond`
"""
bondtype(lat::Lattice, bond::Integer) = lat.bondtypes[bond]

"""
    bonddirection(lat::Lattice, bond::Integer)

return the direction of the `bond` as vector in the Cartesian system
"""
bonddirection(lat::Lattice, bond::Integer) = lat.transvector * lat.bond_dirs[:, bond]

"""
    lattice_bonddirection(lat::Lattice, bond::Integer)

return the direction of the `bond` as vector in the lattice system
"""
lattice_bonddirection(lat::Lattice, bond::Integer) = lat.bond_dirs[:, bond]

"""
    siteL2(lat::Lattice, site::Integer)

return (x+L/2, y+W/2, ...) site
"""
siteL2(lat::Lattice, site::Integer) = lat.site_L2[site]

"""
    siteL4(lat::Lattice, site::Integer)

return (x+L/4, y+W/4, ...) site
"""
siteL4(lat::Lattice, site::Integer) = lat.site_L4[site]

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
    nsites = [1,1]
    nbonds = [1]
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
            nsites, nbonds, sitetypes, bondtypes,
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
    nsites = [L-div(L,2), div(L,2)]
    nbonds = [L-div(L,2), div(L,2)]
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
    Lattice(dim, sz, nsitetypes, nbondtypes,
            nsites, nbonds, sitetypes, bondtypes,
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

    N = L*W
    dim = 2
    sz = [L,W]
    nsitetypes = 2
    nbondtypes = 2
    nsites = zeros(Int, nsitetypes)
    nbonds = zeros(Int, nbondtypes)
    sitetypes = zeros(Int, N)
    bondtypes = zeros(Int, 2N)
    transvector = eye(dim)
    coords = zeros(dim, N)
    bond_dirs = zeros(dim, 2N)
    neighborsites = Vector{Int}[]
    neighborbonds = Vector{Int}[]
    source = zeros(Int,2N)
    target = zeros(Int,2N)
    sitetypes = zeros(Int,N)
    bondtypes = zeros(Int,2N)
    site_L2 = zeros(Int,N)
    site_L4 = zeros(Int,N)
    L2 = div(L,2)
    L4 = div(L,4)
    W2 = div(W,2)
    W4 = div(W,4)
    ns = zeros(Int,4)
    nb = zeros(Int,4)
    @inbounds for s in 1:N
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
        nsites[sitetypes[s]] += 1
        bond_dirs[:, 2s-1] = [1.0, 0.0]
        bond_dirs[:, 2s] = [0.0, 1.0]
        bondtypes[2s-1] = 1 
        bondtypes[2s] = 2 
        nbonds[1] += 1
        nbonds[2] += 1
        site_L2[s] = xy2s(x+L2, y+W2)
        site_L4[s] = xy2s(x+L4, y+W4)
    end
    Lattice(dim, sz, nsitetypes, nbondtypes,
            nsites, nbonds, sitetypes, bondtypes,
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
    N = L*W
    NB = 3N
    nsitetypes = 3
    nbondtypes = 3
    nsites = zeros(Int,nsitetypes)
    nbonds = zeros(Int,nbondtypes)
    sitetypes = zeros(Int,N)
    bondtypes = zeros(Int,NB)
    transvector = [1.0 0.5; 0.0 0.5sqrt(3)]
    coords = zeros(dim, N)
    bond_dirs = zeros(dim, NB)
    neighborsites = Vector{Int}[]
    neighborbonds = Vector{Int}[]
    source = zeros(Int,NB)
    target = zeros(Int,NB)
    site_L2 = zeros(Int, N)
    site_L4 = zeros(Int, N)
    L2 = div(L,2)
    L4 = div(L,4)
    W2 = div(W,2)
    W4 = div(W,4)
    ns = zeros(Int,6)
    nb = zeros(Int,6)
    @inbounds for s in 1:N
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
        nsites[sitetypes[s]] += 1
        bond_dirs[:, 3s-2] .= [1.0,0.0]
        bond_dirs[:, 3s-1] .= [0.0,1.0]
        bond_dirs[:, 3s] .= [-1.0,1.0]
        bondtypes[3s-2] = 1
        bondtypes[3s-1] = 2
        bondtypes[3s] = 3
        nbonds .+= [1,1,1]
        site_L2[s] = xy2s(x+L2, y+W2)
        site_L4[s] = xy2s(x+L4, y+W4)
    end
    Lattice(dim, sz, nsitetypes, nbondtypes,
            nsites, nbonds, sitetypes, bondtypes,
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
    N = L*W*H
    NB = 3N
    nsitetypes = 2
    nbondtypes = 3
    nsites = zeros(Int,nsitetypes)
    nbonds = zeros(Int,nbondtypes)
    sitetypes = zeros(Int,N)
    bondtypes = zeros(Int,NB)
    transvector = eye(dim)
    coords = zeros(dim, N)
    bond_dirs = zeros(dim, NB)
    neighborsites = Vector{Int}[]
    neighborbonds = Vector{Int}[]
    source = zeros(Int,NB)
    target = zeros(Int,NB)
    site_L2 = zeros(Int,N)
    site_L4 = zeros(Int,N)
    L2 = div(L,2)
    L4 = div(L,4)
    W2 = div(W,2)
    W4 = div(W,4)
    H2 = div(H,2)
    H4 = div(H,4)
    ns = zeros(Int,6)
    nb = zeros(Int,6)
    @inbounds for s in 1:N
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
        nsites[sitetypes[s]] += 1
        bondtypes[3s-2] = 1
        bondtypes[3s-1] = 2
        bondtypes[3s] = 3
        bond_dirs[:, 3s-2] = [1.0, 0.0, 0.0]
        bond_dirs[:, 3s-1] = [0.0, 1.0, 0.0]
        bond_dirs[:, 3s-0] = [0.0, 0.0, 1.0]
        nbonds .+= [1,1,1]
        site_L2[s] = xyz2s(x+L2,y+W2,z+H2)
        site_L4[s] = xyz2s(x+L4,y+W4,z+H4)
    end
    Lattice(dim, sz, nsitetypes, nbondtypes,
            nsites, nbonds, sitetypes, bondtypes,
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
