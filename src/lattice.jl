type Lattice
  dim :: Int
  size :: Vector{Int}
  nsites :: Int
  nbonds :: Int
  neighbors :: Matrix{Int}
  source :: Vector{Int}
  target :: Vector{Int}
end

import Base.size

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

return the number of bonds.
"""
numsites(lat::Lattice) = lat.nsites

"""
    numbonds(lat::Lattice)

return the number of bonds.
"""
numbonds(lat::Lattice) = lat.nbonds

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
    chain_lattice(L::Integer)
    
generate chain lattice with length `L`.
"""
function chain_lattice(L::Integer)
  neighbors = zeros(Int,2,L)
  source = zeros(Int,L)
  target = zeros(Int,L)
  for s in 1:L
    neighbors[1,s] = mod1(s+1,L)
    neighbors[2,s] = mod1(s-1,L)
    source[s] = s
    target[s] = neighbors[1,s]
  end
  Lattice(1,[L],L,L,neighbors,source,target)
end

"""
    square_lattice(L::Integer, W::Integer=L)
    
generate square lattice with size `L` \\times `W`.
"""
square_lattice(L::Integer) = square_lattice(L,L)
function square_lattice(L::Integer, W::Integer)
  s2xy(s::Integer) = mod(s-1,L),div(s-1,L)
  xy2s(x::Integer, y::Integer) = mod(y,W)*L + mod(x,L) + 1

  nsites = L*W
  nbonds = 2nsites
  neighbors = zeros(Int,4,nsites)
  source = zeros(Int,nbonds)
  target = zeros(Int,nbonds)
  for s in 1:nsites
    x,y = s2xy(s)
    neighbors[1,s] = xy2s(x+1,y)
    neighbors[2,s] = xy2s(x,y+1)
    neighbors[3,s] = xy2s(x-1,y)
    neighbors[4,s] = xy2s(x,y-1)
    source[2s-1:2s] = s
    target[2s-1:2s] = neighbors[1:2,s]
  end
  Lattice(2,[L,W],nsites,nbonds,neighbors,source,target)
end

"""
    triangular_lattice(L::Integer, W::Integer=L)
    
generate triangular lattice with size `L` \\times `W`.
"""
triangular_lattice(L::Integer) = triangular_lattice(L,L)
function triangular_lattice(L::Integer, W::Integer)
  s2xy(s::Integer) = mod(s-1,L),div(s-1,L)
  xy2s(x::Integer, y::Integer) = mod(y,W)*L + mod(x,L) + 1
  nsites = L*W
  nbonds = 3nsites
  neighbors = zeros(Int,6,nsites)
  source = zeros(Int,nbonds)
  target = zeros(Int,nbonds)
  for s in 1:nsites
    x,y = s2xy(s)
    neighbors[1,s] = xy2s(x+1,y)
    neighbors[2,s] = xy2s(x,y+1)
    neighbors[3,s] = xy2s(x+1,y+1)
    neighbors[4,s] = xy2s(x-1,y)
    neighbors[5,s] = xy2s(x,y-1)
    neighbors[6,s] = xy2s(x-1,y-1)
    source[3s-2:3s] = s
    target[3s-2:3s] = neighbors[1:3,s]
  end
  Lattice(2,[L,W],nsites,nbonds,neighbors,source,target)
end
