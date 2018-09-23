@doc """
    Lattice

# Fields
- `dim :: Int`
    - dimension of lattice
- `size :: Vector{Int}`
    - length of lattice in each dimension
- `nsitetypes :: Int`
    - the number of `sitetype`s
- `nbondtypes :: Int`
    - the number of `bondtype`s
- `sites :: Vector{Vector{Int}}`
    - list of site indecies for each `sitetype`
- `bonds :: Vector{Vector{Int}}`
    - list of bond indecies for each `sitetype`
- `nsites :: Int`
    - the number of sites
- `nbonds :: Int`
    - the number of bonds
- `sitetypes :: Vector{Int}`
    - `sitetype` of each site
- `bondtypes :: Vector{Int}`
    - `bondtype` of each bond
- `transvector :: Matrix{Float64}`
    - lattice vector represented in Cartesian system
    - `@assert size(transvector) == (dim, dim)`
- `site_coords :: Matrix{Float64}`
    - coordinate of each site represented in lattice system
    - `@assert size(site_coords) == (dim, nsites)`
- `bond_dirs :: Matrix{Float64}`
    - displacement of each bond represented in lattice system
    - `@assert size(bond_dirs) == (dim, nbonds)`
- `neighborsites :: Vector{Vector{Int}}`
    - list of indecies of neighbor sites of each site
- `neighborbonds :: Vector{Vector{Int}}`
    - list of indecies of neighbor bonds of each site
- `source :: Vector{Int}`
    - the index of an end site of each bond index
- `target :: Vector{Int}`
    - the index of the other end site of each bond index
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
end

import Base.size
import Distributions.dim

@doc """
    dim(lat::Lattice)
    dim(model::Model)

Returns the dimension of lattice.
"""
dim(lat::Lattice) = lat.dim
dim(model::Model) = dim(model.lat)

@doc """
    size(lat::Lattice, [dim::Integer])
    size(model::Model, [dim::Integer])

Returns the size of lattice.
"""
size(lat::Lattice) = lat.size
size(lat::Lattice, dim::Integer) = lat.size[dim]
size(model::Model) = size(model.lat)
size(model::Model, dim::Integer) = size(model.lat,dim)

@doc """
    sites(lat::Lattice, sitetype::Integer)
    sites(model::Model, sitetype::Integer)

Returns sites with `sitetype`
"""
sites(lat::Lattice, sitetype::Integer) = lat.sites[sitetype]
sites(model::Model, sitetype::Integer) = sites(model.lat, sitetype)

@doc """
    numsites(lat::Lattice)
    numsites(model::Model)

Returns the number of all sites.
"""
numsites(lat::Lattice) = lat.nsites
numsites(model::Model) = numsites(model.lat)

@doc """
    numsites(lat::Lattice, sitetype::Integer)
    numsites(model::Model, sitetype::Integer)

Returns the number of `sitetype` sites.
"""
numsites(lat::Lattice, sitetype::Integer) = length(lat.sites[sitetype])
numsites(model::Model, sitetype::Integer) = numsites(model.lat, sitetype)

@doc """
    bonds(lat::Lattice, bondtype::Integer)
    bonds(model::Model, bondtype::Integer)

Returns bonds with `bondtype`
"""
bonds(lat::Lattice, bondtype::Integer) = lat.bonds[bondtype]
bonds(model::Model, bondtype::Integer) = bonds(model.lat, bondtype)

@doc """
    numbonds(lat::Lattice)
    numbonds(model::Model)

Returns the number of all bonds.
"""
numbonds(lat::Lattice) = lat.nbonds
numbonds(model::Model) = numbonds(model.lat)

@doc """
    numbonds(lat::Lattice, bondtype::Integer)
    numbonds(model::Model, bondtype::Integer)

Returns the number of `bondtype` bonds.
"""
numbonds(lat::Lattice, bondtype::Integer) = length(lat.bonds[bondtype])
numbonds(model::Model, bondtype::Integer) = numbonds(model.lat, bondtype)

@doc """
    numsitetypes(lat::Lattice)
    numsitetypes(model::Model)

Returns the number of sitetypes.
"""
numsitetypes(lat::Lattice) = lat.nsitetypes
numsitetypes(model::Model) = numsitetypes(model.lat)

@doc """
    numbondtypes(lat::Lattice)
    numbondtypes(model::Model)

Returns the number of bondtypes.
"""
numbondtypes(lat::Lattice) = lat.nbondtypes
numbondtypes(model::Model) = numbondtypes(model.lat)

@doc """
    neighborsites(lat::Lattice, site::Integer)
    neighborsites(model::Model, site::Integer)

Returns the neighbor sites of `site`.
"""
neighborsites(lat::Lattice, site::Integer) = lat.neighborsites[site]
neighborsites(model::Model, site::Integer) = neighborsites(model.lat, site)

@doc """
    neighborbonds(lat::Lattice, site::Integer)
    neighborbonds(model::Model, site::Integer)

Returns the neighbor bonds of `site`.
"""
neighborbonds(lat::Lattice, site::Integer) = lat.neighborbonds[site]
neighborbonds(model::Model, site::Integer) = neighborbonds(model.lat, site)

@doc """
    neighbors(lat::Lattice, site::Integer)
    neighbors(model::Model, site::Integer)

Returns the neighbor sites and bonds of `site`.
"""
neighbors(lat::Lattice, site::Integer) = zip(neighborsites(lat,site), neighborbonds(lat,site))
neighbors(model::Model, site::Integer) = neighbors(model.lat, site)

@doc """
    source(lat::Lattice, bond::Integer)
    source(model::Model, bond::Integer)

Returns the source site of `bond`.
"""
source(lat::Lattice, bond::Integer) = lat.source[bond]
source(model::Model, bond::Integer) = source(model.lat, bond)

@doc """
    target(lat::Lattice, bond::Integer)
    target(model::Model, bond::Integer)

Returns the target site of `bond`.
"""
target(lat::Lattice, bond::Integer) = lat.target[bond]
target(model::Model, bond::Integer) = target(model.lat, bond)

@doc """
    sitecoordinate(lat::Lattice, site::Integer)
    sitecoordinate(model::Model, site::Integer)

Returns the coordinate of the `site` in the Cartesian system
"""
sitecoordinate(lat::Lattice, site::Integer) = lat.transvector * lat.site_coords[:,site]
sitecoordinate(model::Model, site::Integer) = sitecoordinate(model.lat, site)

@doc """
    lattice_sitecoordinate(lat::Lattice, site::Integer)
    lattice_sitecoordinate(model::Model, site::Integer)

Returns the coordinate of the `site` in the lattice system
"""
lattice_sitecoordinate(lat::Lattice, site::Integer) = lat.site_coords[:, site]
lattice_sitecoordinate(model::Model, site::Integer) = lattice_sitecoordinate(model.lat, site)

@doc """
    sitetype(lat::Lattice, site::Integer)
    sitetype(model::Model, site::Integer)

Returns the type of `site`
"""
sitetype(lat::Lattice, site::Integer) = lat.sitetypes[site]
sitetype(model::Model, site::Integer) = sitetype(model.lat, site)

@doc """
    bondtype(lat::Lattice, bond::Integer)
    bondtype(model::Model, bond::Integer)

Returns the type of `bond`
"""
bondtype(lat::Lattice, bond::Integer) = lat.bondtypes[bond]
bondtype(model::Model, bond::Integer) = bondtype(model.lat, bond)

@doc """
    bonddirection(lat::Lattice, bond::Integer)
    bonddirection(model::Model, bond::Integer)

Returns the direction of the `bond` as vector in the Cartesian system
"""
bonddirection(lat::Lattice, bond::Integer) = lat.transvector * lat.bond_dirs[:, bond]
bonddirection(model::Model, bond::Integer) = bonddirection(model.lat, bond)

@doc """
    lattice_bonddirection(lat::Lattice, bond::Integer)
    lattice_bonddirection(model::Model, bond::Integer)

Returns the direction of the `bond` as vector in the lattice system
"""
lattice_bonddirection(lat::Lattice, bond::Integer) = lat.bond_dirs[:, bond]
lattice_bonddirection(model::Model, bond::Integer) = lattice_bonddirection(model.lat, bond)

