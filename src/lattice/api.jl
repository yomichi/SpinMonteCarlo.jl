import Base.size
import Distributions.dim
import LinearAlgebra.norm

# struct Model
#     lat :: Lattice
# end

@doc """
    dim(lat::Lattice)
    dim(model::Model)

Returns the dimension of lattice.
"""
@inline dim(lat::Lattice) = size(lat.latvec,1)
@inline dim(model::Model) = dim(model.lat)

@doc """
    size(lat::Lattice, [dim::Integer])
    size(model::Model, [dim::Integer])

Returns the size of lattice.
"""
@inline size(lat::Lattice) = lat.L
@inline size(lat::Lattice, dim::Integer) = size(lat)[dim]
@inline size(model::Model) = size(model.lat)
@inline size(model::Model, dim::Integer) = size(model.lat,dim)

@doc """
    sites(lat::Lattice, [sitetype::Integer])
    sites(model::Model, [sitetype::Integer])

Returns an iterator over sites with `sitetype` (if omitted, over all sites)
"""
@inline sites(lat::Lattice) = lat.sites
@inline sites(model::Model) = sites(model.lat)
@inline sites(lat::Lattice, sitetype::Integer) = lat.siteswithtype[sitetype]
@inline sites(model::Model, sitetype::Integer) = sites(model.lat, sitetype)

@doc """
    numsites(lat::Lattice)
    numsites(model::Model)

Returns the number of all sites.
"""
@inline numsites(lat::Lattice) = length(lat.sites)
@inline numsites(model::Model) = numsites(model.lat)

@doc """
    numsites(lat::Lattice, sitetype::Integer)
    numsites(model::Model, sitetype::Integer)

Returns the number of `sitetype` sites.
"""
@inline numsites(lat::Lattice, sitetype::Integer) = length(lat.siteswithtype[sitetype])
@inline numsites(model::Model, sitetype::Integer) = numsites(model.lat, sitetype)

@doc """
    bonds(lat::Lattice, [bondtype::Integer])
    bonds(model::Model, [bondtype::Integer])

Returns an iterator over bonds with `bondtype` (if omitted, over all bonds)
"""
@inline bonds(lat::Lattice) = lat.bonds
@inline bonds(model::Model) = bonds(model.lat)
@inline bonds(lat::Lattice, bondtype::Integer) = lat.bondswithtype[bondtype]
@inline bonds(model::Model, bondtype::Integer) = bonds(model.lat, bondtype)

@doc """
    numbonds(lat::Lattice)
    numbonds(model::Model)

Returns the number of all bonds.
"""
@inline numbonds(lat::Lattice) = length(lat.bonds)
@inline numbonds(model::Model) = numbonds(model.lat)

@doc """
    numbonds(lat::Lattice, bondtype::Integer)
    numbonds(model::Model, bondtype::Integer)

Returns the number of `bondtype` bonds.
"""
@inline numbonds(lat::Lattice, bondtype::Integer) = length(lat.bondswithtype[bondtype])
@inline numbonds(model::Model, bondtype::Integer) = numbonds(model.lat, bondtype)

@doc """
    numsitetypes(lat::Lattice)
    numsitetypes(model::Model)

Returns the number of sitetypes.
"""
@inline numsitetypes(lat::Lattice) = length(lat.siteswithtype)
@inline numsitetypes(model::Model) = numsitetypes(model.lat)

@doc """
    numbondtypes(lat::Lattice)
    numbondtypes(model::Model)

Returns the number of bondtypes.
"""
@inline numbondtypes(lat::Lattice) = length(lat.bondswithtype)
@inline numbondtypes(model::Model) = numbondtypes(model.lat)

@doc """
    neighborsites(lat::Lattice, site::Integer)
    neighborsites(model::Model, site::Integer)

Returns the neighbor sites of `site`.
"""
@inline neighborsites(lat::Lattice, site::Integer) = neighbors(lat, site)
@inline neighborsites(model::Model, site::Integer) = neighborsites(model.lat, site)

@doc """
    source(lat::Lattice, bond::Integer)
    source(model::Model, bond::Integer)

Returns the source site of `bond`.
"""
@inline source(lat::Lattice, bond::Integer) = lat.bonds[bond].source
@inline source(model::Model, bond::Integer) = source(model.lat, bond)

@doc """
    target(lat::Lattice, bond::Integer)
    target(model::Model, bond::Integer)

Returns the target site of `bond`.
"""
@inline target(lat::Lattice, bond::Integer) = lat.bonds[bond].target
@inline target(model::Model, bond::Integer) = target(model.lat, bond)

@doc """
    sitecoordinate(lat::Lattice, site::Integer)
    sitecoordinate(model::Model, site::Integer)

Returns the coordinate of the `site` in the Cartesian system
"""
@inline sitecoordinate(lat::Lattice, site::Integer) = lat.sites[site].coord
@inline sitecoordinate(model::Model, site::Integer) = sitecoordinate(model.lat, site)

@doc """
    cellcoordinate(lat::Lattice, site::Integer)
    cellcoordinate(model::Model, site::Integer)

Returns the coordinate of the `cell` including `site` in the Lattice system
"""
@inline lattice_sitecoordinate(lat::Lattice, site::Integer) = lat.sites[site].cellcoord
@inline lattice_sitecoordinate(model::Model, site::Integer) = lattice_sitecoordinate(model.lat, site)

@doc """
    sitetype(lat::Lattice, site::Integer)
    sitetype(model::Model, site::Integer)

Returns the type of `site`
"""
@inline sitetype(lat::Lattice, site::Integer) = lat.sites[site].sitetype
@inline sitetype(model::Model, site::Integer) = sitetype(model.lat, site)

@doc """
    bondtype(lat::Lattice, bond::Integer)
    bondtype(model::Model, bond::Integer)

Returns the type of `bond`.
"""
@inline bondtype(lat::Lattice, bond::Integer) = lat.bonds[bond].bondtype
@inline bondtype(model::Model, bond::Integer) = bondtype(model.lat, bond)

@doc """
    bonddirection(lat::Lattice, bond::Integer)
    bonddirection(model::Model, bond::Integer)

Returns the unnormalized direction of the `bond` as vector in the Cartesian system
"""
bonddirection(lat::Lattice, bond::Integer) = lat.bonds[bond].direction
bonddirection(model::Model, bond::Integer) = bonddirection(model.lat, bond)
