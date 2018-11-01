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
@inline dim(lat::Lattice) = get_prop(lat,:dim)
@inline dim(model::Model) = dim(model.lat)

@doc """
    size(lat::Lattice, [dim::Integer])
    size(model::Model, [dim::Integer])

Returns the size of lattice.
"""
@inline size(lat::Lattice) = get_prop(lat, :size)
@inline size(lat::Lattice, dim::Integer) = size(lat)[dim]
@inline size(model::Model) = size(model.lat)
@inline size(model::Model, dim::Integer) = size(model.lat,dim)

@doc """
    sites(lat::Lattice, [sitetype::Integer])
    sites(model::Model, [sitetype::Integer])

Returns an iterator over sites with `sitetype` (if omitted, over all sites)
"""
@inline sites(lat::Lattice) = vertices(lat)
@inline sites(model::Model) = sites(model.lat)
@inline sites(lat::Lattice, sitetype::Integer) = (site for site in vertices(lat) if get_prop(lat, site, :sitetype)==sitetype)
@inline sites(model::Model, sitetype::Integer) = sites(model.lat, sitetype)

@doc """
    numsites(lat::Lattice)
    numsites(model::Model)

Returns the number of all sites.
"""
@inline numsites(lat::Lattice) = get_prop(lat, :numsites)
@inline numsites(model::Model) = numsites(model.lat)

@doc """
    numsites(lat::Lattice, sitetype::Integer)
    numsites(model::Model, sitetype::Integer)

Returns the number of `sitetype` sites.
"""
@inline numsites(lat::Lattice, sitetype::Integer) = get_prop(lat, :nss)[sitetype]
@inline numsites(model::Model, sitetype::Integer) = numsites(model.lat, sitetype)

@doc """
    bonds(lat::Lattice, [bondtype::Integer])
    bonds(model::Model, [bondtype::Integer])

Returns an iterator over bonds with `bondtype` (if omitted, over all bonds)
"""
@inline bonds(lat::Lattice) = edges(lat)
@inline bonds(model::Model) = bonds(model.lat)
@inline bonds(lat::Lattice, bondtype::Integer) = (bond for bond in edges(lat) if get_prop(lat, bond, :bondtype)==bondtype)
@inline bonds(model::Model, bondtype::Integer) = bonds(model.lat, bondtype)

@doc """
    numbonds(lat::Lattice)
    numbonds(model::Model)

Returns the number of all bonds.
"""
@inline numbonds(lat::Lattice) = get_prop(lat, :numbonds)
@inline numbonds(model::Model) = numbonds(model.lat)

@doc """
    numbonds(lat::Lattice, bondtype::Integer)
    numbonds(model::Model, bondtype::Integer)

Returns the number of `bondtype` bonds.
"""
@inline numbonds(lat::Lattice, bondtype::Integer) = get_prop(lat, :nbs)[bondtype]
@inline numbonds(model::Model, bondtype::Integer) = numbonds(model.lat, bondtype)

@doc """
    numsitetypes(lat::Lattice)
    numsitetypes(model::Model)

Returns the number of sitetypes.
"""
@inline numsitetypes(lat::Lattice) = get_prop(lat, :numsitetypes)
@inline numsitetypes(model::Model) = numsitetypes(model.lat)

@doc """
    numbondtypes(lat::Lattice)
    numbondtypes(model::Model)

Returns the number of bondtypes.
"""
@inline numbondtypes(lat::Lattice) = get_prop(lat, :numbondtypes)
@inline numbondtypes(model::Model) = numbondtypes(model.lat)

@doc """
    neighborsites(lat::Lattice, site::Integer)
    neighborsites(model::Model, site::Integer)

Returns the neighbor sites of `site`.
"""
@inline neighborsites(lat::Lattice, site::Integer) = neighbors(lat, site)
@inline neighborsites(model::Model, site::Integer) = neighborsites(model.lat, site)

@doc """
    source(lat::Lattice, bond::Edge)
    source(model::Model, bond::Edge)
    source(lat::Lattice, s1::Integer, s2::Integer)
    source(model::Model, s1::Integer, s2::Integer)

Returns the source site of `bond` or `Edge(s1,s2)`.
"""
@inline source(lat::Lattice, bond::Edge) = get_prop(lat, bond, :source)
@inline source(model::Model, bond::Edge) = source(model.lat, bond)
@inline source(lat::Lattice, s1::Integer, s2::Integer) = source(lat, Edge(s1,s2))
@inline source(model::Model, s1::Integer, s2::Integer) = source(model.lat, Edge(s1,s2))

@doc """
    target(lat::Lattice, bond::Edge)
    target(model::Model, bond::Edge)
    target(lat::Lattice, s1::Integer, s2::Integer)
    target(model::Model, s1::Integer, s2::Integer)

Returns the target site of `bond` or `Edge(s1,s2).
"""
@inline target(lat::Lattice, bond::Edge) = get_prop(lat, bond, :target)
@inline target(model::Model, bond::Edge) = target(model.lat, bond)
@inline target(lat::Lattice, s1::Integer, s2::Integer) = target(lat, Edge(s1,s2))
@inline target(model::Model, s1::Integer, s2::Integer) = target(model.lat, Edge(s1,s2))

@doc """
    sitecoordinate(lat::Lattice, site::Integer)
    sitecoordinate(model::Model, site::Integer)

Returns the coordinate of the `site` in the Cartesian system
"""
@inline sitecoordinate(lat::Lattice, site::Integer) = get_prop(lat, site, :coord)
@inline sitecoordinate(model::Model, site::Integer) = sitecoordinate(model.lat, site)

@doc """
    cellcoordinate(lat::Lattice, site::Integer)
    cellcoordinate(model::Model, site::Integer)

Returns the coordinate of the `cell` including `site` in the Lattice system
"""
@inline lattice_sitecoordinate(lat::Lattice, site::Integer) = get_prop(lat, site, :cellcoord)
@inline lattice_sitecoordinate(model::Model, site::Integer) = lattice_sitecoordinate(model.lat, site)

@doc """
    sitetype(lat::Lattice, site::Integer)
    sitetype(model::Model, site::Integer)

Returns the type of `site`
"""
@inline sitetype(lat::Lattice, site::Integer) = get_prop(lat, site, :sitetype)
@inline sitetype(model::Model, site::Integer) = sitetype(model.lat, site)

@doc """
    bondtype(lat::Lattice, bond::Edge)
    bondtype(model::Model, bond::Edge)
    bondtype(lat::Lattice, s1::Integer, s2::Integer)
    bondtype(model::Model, s1::Integer, s2::Integer)

Returns the type of `bond` or `Edge(s1,s2)`
"""
@inline bondtype(lat::Lattice, bond::Edge) = get_prop(lat, bond, :bondtype)
@inline bondtype(model::Model, bond::Edge) = bondtype(model.lat, bond)
@inline bondtype(lat::Lattice, s1::Integer, s2::Integer) = bondtype(lat, Edge(s1,s2))
@inline bondtype(model::Model, s1::Integer, s2::Integer) = bondtype(model.lat, Edge(s2,s2))

@doc """
    bonddirection(lat::Lattice, bond::Edge)
    bonddirection(model::Model, bond::Edge)
    bonddirection(lat::Lattice, s1::Integer, s2::Integer)
    bonddirection(model::Model, s1::Integer, s2::Integer)

Returns the unnormalized direction of the `bond` as vector in the Cartesian system
"""
bonddirection(lat::Lattice, bond::Edge) = get_prop(lat, bond, :bonddirection)
bonddirection(model::Model, bond::Edge) = bonddirection(model.lat, bond)
bonddirection(lat::Lattice, s1::Integer, s2::Integer) = bonddirection(lat, Edge(s1,s2))
bonddirection(model::Model, s1::Integer, s2::Integer) = bonddirection(model.lat, Edge(s1,s2))
