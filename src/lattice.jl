"""
    Lattice
"""
mutable struct Lattice
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

