import Base.convert

mutable struct Site
    id :: Int
    sitetype :: Int
    neighborsites :: Vector{Int}
    neighborbonds :: Vector{Int}
    coord :: Vector{Float64}
    localsite :: Int
    cellcoord :: Vector{Int}
end
convert(::Type{Int}, s::Site) = s.id

mutable struct Bond
    id :: Int
    bondtype :: Int
    source :: Int
    target :: Int
    direction :: Vector{Float64}
end
convert(::Type{Int}, b::Bond) = b.id

mutable struct Lattice
    latticevector :: Matrix{Float64}
    L :: Vector{Int}
    sites :: Vector{Site}
    siteswithtype :: Vector{SubArray{Site,1,Vector{Site},Tuple{Vector{Int}},false}}
    bonds :: Vector{Bond}
    bondswithtype :: Vector{SubArray{Bond,1,Vector{Bond},Tuple{Vector{Int}},false}}

    function Lattice(latvec,L,sites,bonds)
        nsitetypes = 0
        swt = Vector{Int}[]
        for s in sites
            st = s.sitetype
            while st > nsitetypes
                push!(swt, Int[])
                nsitetypes += 1
            end
            push!(swt[st], s.id)
        end
        siteswithtype = [view(sites, swt[st]) for st in 1:nsitetypes]

        nbondtypes = 0
        bwt = Vector{Int}[]
        for b in bonds
            bt = b.bondtype
            while bt > nbondtypes
                push!(bwt, Int[])
                nbondtypes += 1
            end
            push!(bwt[bt], b.id)

            push!(sites[source(b)].neighborsites, target(b))
            push!(sites[source(b)].neighborbonds, b.id)
            push!(sites[target(b)].neighborsites, source(b))
            push!(sites[target(b)].neighborbonds, b.id)
        end
        bondswithtype = [view(bonds, bwt[bt]) for bt in 1:nbondtypes]

        return new(latvec, L,
                   sites, siteswithtype,
                   bonds, bondswithtype,
                  )
    end
end

include("parameter.jl")
include("api.jl")
include("generate.jl")
include("standard.jl")

export LatticeParameter
export generatelattice
