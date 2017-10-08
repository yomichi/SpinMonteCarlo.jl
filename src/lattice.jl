type Lattice
    dim :: Int
    size :: Vector{Int}
    nsites :: Int
    nbonds :: Int
    nsitetypes :: Int
    nbondtypes :: Int
    site_coords :: Matrix{Float64}
    bond_dirs :: Matrix{Float64}
    neighbors :: Matrix{Int}
    source :: Vector{Int}
    target :: Vector{Int}
    sitetypes :: Vector{Int}
    bondtypes :: Vector{Int}
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

return the number of sites.
"""
numsites(lat::Lattice) = lat.nsites

"""
    numbonds(lat::Lattice)

return the number of bonds.
"""
numbonds(lat::Lattice) = lat.nbonds

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
    neighbors(lat::Lattice, site::Integer)

return the array of neighbor sites of `site`.
"""
neighbors(lat::Lattice, site::Integer) = lat.neighbors[:,site]

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
    coordinate(lat::Lattice, site::Integer)

return the coordinate of the `site`
"""
coordinate(lat::Lattice, site::Integer) = lat.site_coords[:, site]

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
    bonddirectory(lat::Lattice, bond::Integer)

return the directory of the `bond` as vector
"""
bonddirectory(lat::Lattice, bond::Integer) = lat.bond_dirs[:, bond]

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
    coords = zeros(dim, 2)
    bond_dirs = zeros(dim, 1)
    neighbors = zeros(Int,1,2)
    source = zeros(Int,1)
    target = zeros(Int,1)
    site_L2 = zeros(Int,2)
    site_L4 = zeros(Int,2)

    coords[1,1] = 0.0
    coords[1,2] = 1.0
    bond_dirs[1,1] = 1.0
    neighbors[1,1] = 2
    neighbors[1,2] = 1
    source[1] = 1
    target[1] = 2
    site_L2[1] = 2
    site_L2[2] = 1
    site_L4[1] = 2
    site_L4[2] = 1

    Lattice(dim,[2],2,1,coords, bond_dirs, neighbors,source,target, site_L2, site_L4)
end
dimer_lattice(params::Dict) = dimer_lattice()

"""
    chain_lattice(L::Integer)
    chain_lattice(params::Dict)
    
generate chain lattice with length `L` or params["L"].
"""
function chain_lattice(L::Integer)
    dim = 1
    coords = zeros(dim, L)
    bond_dirs = zeros(dim, L)
    neighbors = zeros(Int,2,L)
    source = zeros(Int,L)
    target = zeros(Int,L)
    sitetypes = zeros(Int,L)
    bondtypes = zeros(Int,L)
    site_L2 = zeros(Int,L)
    site_L4 = zeros(Int,L)
    for s in 1:L
        coords[1,s] = (s-1)
        bond_dirs[1,s] = 1.0
        neighbors[1,s] = mod1(s+1,L)
        neighbors[2,s] = mod1(s-1,L)
        source[s] = s
        target[s] = neighbors[1,s]
        sitetypes[s] = ifelse(isodd(s),1,2)
        bondtypes[s] = ifelse(isodd(s),1,2)
        site_L2[s] = mod1(s+div(L,2),L)
        site_L4[s] = mod1(s+div(L,4),L)
    end
    Lattice(dim, [L], L, L, 2, 2, coords, bond_dirs,
            neighbors, source, target,
            sitetypes, bondtypes, site_L2, site_L4)
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
    nsites = L*W
    nbonds = 2nsites
    coords = zeros(dim, nsites)
    bond_dirs = zeros(dim, nbonds)
    neighbors = zeros(Int,4,nsites)
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
    @inbounds for s in 1:nsites
        x,y = s2xy(s)
        neighbors[1,s] = xy2s(x+1,y)
        neighbors[2,s] = xy2s(x,y+1)
        neighbors[3,s] = xy2s(x-1,y)
        neighbors[4,s] = xy2s(x,y-1)
        source[2s-1:2s] = s
        target[2s-1:2s] = neighbors[1:2,s]
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
    Lattice(dim, [L,W], nsites, nbonds, 2, 2, coords, bond_dirs,
            neighbors, source, target,
            sitetypes, bondtypes, site_L2, site_L4)
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
    nsites = L*W
    nbonds = 3nsites
    coords = zeros(dim, nsites)
    bond_dirs = zeros(dim, nbonds)
    neighbors = zeros(Int,6,nsites)
    source = zeros(Int,nbonds)
    target = zeros(Int,nbonds)
    sitetypes = zeros(Int,nsites)
    bondtypes = zeros(Int,nbonds)
    site_L2 = zeros(Int, nsites)
    site_L4 = zeros(Int, nsites)
    L2 = div(L,2)
    L4 = div(L,4)
    W2 = div(W,2)
    W4 = div(W,4)
    @inbounds for s in 1:nsites
        x,y = s2xy(s)
        neighbors[1,s] = xy2s(x+1,y)
        neighbors[2,s] = xy2s(x,y+1)
        neighbors[3,s] = xy2s(x+1,y+1)
        neighbors[4,s] = xy2s(x-1,y)
        neighbors[5,s] = xy2s(x,y-1)
        neighbors[6,s] = xy2s(x-1,y-1)
        source[3s-2:3s] = s
        target[3s-2:3s] = neighbors[1:3,s]
        coords[:,s] = x*ex + y*ey
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
        bond_dirs[:, 3s-2] = ex
        bond_dirs[:, 3s-1] = ey
        bond_dirs[:, 3s] = ex+ey
        bondtypes[3s-2] = 1
        bondtypes[3s-1] = 2
        bondtypes[3s] = 3
        site_L2[s] = xy2s(x+L2, y+W2)
        site_L4[s] = xy2s(x+L4, y+W4)
    end
    Lattice(dim, [L,W], nsites, nbonds, 3, 3,
            coords, bond_dirs, neighbors, source, target,
            sitetypes, bondtypes, site_L2,site_L4)
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
    nsites = L*W*H
    nbonds = 3nsites
    coords = zeros(dim, nsites)
    bond_dirs = zeros(dim, nbonds)
    neighbors = zeros(Int,6,nsites)
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
    H2 = div(H,2)
    H4 = div(H,4)
    @inbounds for s in 1:nsites
        x,y,z = s2xyz(s)
        neighbors[1,s] = xyz2s(x+1,y,z)
        neighbors[2,s] = xyz2s(x,y+1,z)
        neighbors[3,s] = xyz2s(x,y,z+1)
        neighbors[4,s] = xyz2s(x-1,y,z)
        neighbors[5,s] = xyz2s(x,y-1,z)
        neighbors[6,s] = xyz2s(x,y,z-1)
        source[3s-2:3s] = s
        target[3s-2:3s] = neighbors[1:3,s]
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
        site_L4[s] = xyz4s(x+L4,y+W4,z+H4)
    end
    Lattice(dim, [L,W,H], nsites, nbonds, 2, 3,
            coords, bond_dirs, neighbors, source, target,
            sitetypes, bondtypes, site_L2, site_L4)
end
function cubic_lattice(params::Dict)
    L = params["L"]
    W = get(params, "W", L)
    H = get(params, "H", W)
    return cubic_lattice(L, W, H)
end
