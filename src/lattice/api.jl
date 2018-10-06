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
dim(lat::Lattice) = get_prop(lat,:dim)
dim(model::Model) = dim(model.lat)

@doc """
    size(lat::Lattice, [dim::Integer])
    size(model::Model, [dim::Integer])

Returns the size of lattice.
"""
size(lat::Lattice) = get_prop(lat, :size)
size(lat::Lattice, dim::Integer) = size(lat)[dim]
size(model::Model) = size(model.lat)
size(model::Model, dim::Integer) = size(model.lat,dim)

@doc """
    sites(lat::Lattice, sitetype::Integer)
    sites(model::Model, sitetype::Integer)

Returns an iterator over sites with `sitetype`
"""
sites(lat::Lattice, sitetype::Integer) = (site for site in vertices(lat) if get_prop(lat, site, :sitetype)==sitetype)
sites(model::Model, sitetype::Integer) = sites(model.lat, sitetype)

@doc """
    numsites(lat::Lattice)
    numsites(model::Model)

Returns the number of all sites.
"""
numsites(lat::Lattice) = get_prop(lat, :numsites)
numsites(model::Model) = numsites(model.lat)

@doc """
    numsites(lat::Lattice, sitetype::Integer)
    numsites(model::Model, sitetype::Integer)

Returns the number of `sitetype` sites.
"""
numsites(lat::Lattice, sitetype::Integer) = get_prop(g, :nss)[sitetype]
numsites(model::Model, sitetype::Integer) = numsites(model.lat, sitetype)

@doc """
    bonds(lat::Lattice, bondtype::Integer)
    bonds(model::Model, bondtype::Integer)

Returns bonds with `bondtype`
"""
bonds(lat::Lattice, bondtype::Integer) = (bond for bond in edges(lat) if get_prop(lat, bond, :bondtype)==bondtype)
bonds(model::Model, bondtype::Integer) = bonds(model.lat, bondtype)

@doc """
    numbonds(lat::Lattice)
    numbonds(model::Model)

Returns the number of all bonds.
"""
numbonds(lat::Lattice) = get_prop(lat, :numbonds)
numbonds(model::Model) = numbonds(model.lat)

@doc """
    numbonds(lat::Lattice, bondtype::Integer)
    numbonds(model::Model, bondtype::Integer)

Returns the number of `bondtype` bonds.
"""
numbonds(lat::Lattice, bondtype::Integer) = get_prop(lat, :nbs)[bondtype]
numbonds(model::Model, bondtype::Integer) = numbonds(model.lat, bondtype)

@doc """
    numsitetypes(lat::Lattice)
    numsitetypes(model::Model)

Returns the number of sitetypes.
"""
numsitetypes(lat::Lattice) = get_prop(lat, :numsitetypes)
numsitetypes(model::Model) = numsitetypes(model.lat)

@doc """
    numbondtypes(lat::Lattice)
    numbondtypes(model::Model)

Returns the number of bondtypes.
"""
numbondtypes(lat::Lattice) = get_prop(lat, :numbondtypes)
numbondtypes(model::Model) = numbondtypes(model.lat)

@doc """
    neighborsites(lat::Lattice, site::Integer)
    neighborsites(model::Model, site::Integer)

Returns the neighbor sites of `site`.
"""
neighborsites(lat::Lattice, site::Integer) = neighbors(lat, site)
neighborsites(model::Model, site::Integer) = neighborsites(model.lat, site)

@doc """
    source(lat::Lattice, bond::Edge)
    source(model::Model, bond::Edge)

Returns the source site of `bond`.
"""
source(lat::Lattice, bond::Edge) = get_prop(lat, bond, :source)
source(model::Model, bond::Edge) = source(model.lat, bond)

@doc """
    target(lat::Lattice, bond::Edge)
    target(model::Model, bond::Edge)

Returns the target site of `bond`.
"""
target(lat::Lattice, bond::Edge) = get_prop(lat, bond, :target)
target(model::Model, bond::Edge) = target(model.lat, bond)

@doc """
    sitecoordinate(lat::Lattice, site::Integer)
    sitecoordinate(model::Model, site::Integer)

Returns the coordinate of the `site` in the Cartesian system
"""
sitecoordinate(lat::Lattice, site::Integer) = get_prop(lat, site, :coord)
sitecoordinate(model::Model, site::Integer) = sitecoordinate(model.lat, site)

@doc """
    cellcoordinate(lat::Lattice, site::Integer)
    cellcoordinate(model::Model, site::Integer)

Returns the coordinate of the `cell` including `site` in the Lattice system
"""
lattice_sitecoordinate(lat::Lattice, site::Integer) = get_prop(lat, site, :cellcoord)
lattice_sitecoordinate(model::Model, site::Integer) = lattice_sitecoordinate(model.lat, site)

@doc """
    sitetype(lat::Lattice, site::Integer)
    sitetype(model::Model, site::Integer)

Returns the type of `site`
"""
sitetype(lat::Lattice, site::Integer) = get_prop(lat, site, :sitetype)
sitetype(model::Model, site::Integer) = sitetype(model.lat, site)

@doc """
    bondtype(lat::Lattice, bond::Edge)
    bondtype(model::Model, bond::Edge)

Returns the type of `bond`
"""
bondtype(lat::Lattice, bond::Edge) = get_prop(lat, bond, :bondtype)
bondtype(model::Model, bond::Edge) = bondtype(model.lat, bond)

@doc """
    bonddirection(lat::Lattice, bond::Edge)
    bonddirection(model::Model, bond::Edge)

Returns the direction of the `bond` as vector in the Cartesian system
"""
bonddirection(lat::Lattice, bond::Edge) = get_prop(lat, bond, :bonddirection)
bonddirection(model::Model, bond::Edge) = bonddirection(model.lat, bond)
