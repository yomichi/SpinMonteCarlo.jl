mutable struct Site
    id :: Int
    sitetype :: Int
    neighborsites :: Vector{Int}
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
end
Bond() = Bond(0,0,0,0)

mutable struct Lattice
    latticevector :: Matrix{Float64}
    L :: Vector{Int}
    sites :: Vector{Site}
    siteswithtype :: Vector{SubArray{Site,1,Vector{Site},Tuple{Vector{Site}},false}}
    bonds :: Vector{Bond}
    bondswithtype :: Vector{SubArray{Bond,1,Vector{Bond},Tuple{Vector{Bond}},false}}

    function Lattice(latvec,L,sites,bonds)
        nsitetypes = 1
        for (i,s) in enumerate(sites)
            s.id = i
            st = s.sitetype
            swt = Vector{Int}[]
            while st > nsitetypes
                push!(swt, Int[])
                nsitetypes += 1
            end
            push!(swt[st], i)
        end
        siteswithtype = [view(sites, swt[st]) for st in 1:nsitetypes]

        nbondtypes = 1
        for (i,b) in enumerate(bonds)
            b.id = i
            bt = b.bondtype
            bwt = Vector{Int}[]
            while bt > nbondtypes
                push!(bwt, Int[])
                nbondtypes += 1
            end
            push!(bwt[bt], i)
        end
        bondswithtype = [view(bonds, bwt[bt]) for st in 1:nbondtypes]

        return Lattice(latvec, L,
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
