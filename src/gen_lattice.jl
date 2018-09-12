@doc doc"""
    dimer_lattice()
    dimer_lattice(param::Dict)
    
Generates dimer lattice (indeed, this is not a "lattice")
"""
function dimer_lattice()
    dim = 1
    sz = [2]
    nsitetypes = 2
    nbondtypes = 1
    sites = Vector{Int}[zeros(Int,0) for i in 1:nsitetypes]
    bonds = Vector{Int}[zeros(Int,0) for i in 1:nbondtypes]
    nsites = 2
    nbonds = 1
    sitetypes = [1,2]
    bondtypes = [1]
    transvector = @compat Matrix(1.0I, dim, dim)
    coords = zeros(dim, 2)
    bond_dirs = zeros(dim, 1)
    neighborsites = [[2],[1]]
    neighborbonds = [[1],[1]]
    source = zeros(Int,1)
    target = zeros(Int,1)

    sites[1] = [1]
    sites[2] = [2]
    bonds[1] = [1]
    coords[1,1] = 0.0
    coords[1,2] = 1.0
    bond_dirs[1,1] = 1.0
    source[1] = 1
    target[1] = 2

    Lattice(dim, sz, nsitetypes, nbondtypes,
            sites, bonds, nsites, nbonds, sitetypes, bondtypes,
            transvector, coords, bond_dirs,
            neighborsites, neighborbonds,
            source,target,
           )
end
dimer_lattice(param::Dict) = dimer_lattice()

@doc doc"""
    chain_lattice(L::Integer)
    chain_lattice(param::Dict)
    
Generates chain lattice with length `L = param["L"]`

`nsitetypes` is 2: all `1` (`2`) sites do not connect to another `1` (`2`) site.

`nbondtypes` is 2: `1` bond connects `2n-1` and `2n` sites and `2` bond connects `2n` and `2n+1` sites.
"""
function chain_lattice(L::Integer)
    dim = 1
    sz = [L]
    nsitetypes = 2
    nbondtypes = 2
    nsites = L
    nbonds = L
    sitetypes = zeros(Int,L)
    bondtypes = zeros(Int,L)
    transvector = @compat Matrix(1.0I, dim, dim)
    coords = zeros(dim, L)
    bond_dirs = zeros(dim, L)
    neighborsites = Vector{Int}[]
    neighborbonds = Vector{Int}[]
    source = zeros(Int,L)
    target = zeros(Int,L)

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
    end
    sites = Vector{Int}[zeros(Int,0) for i in 1:nsitetypes]
    bonds = Vector{Int}[zeros(Int,0) for i in 1:nbondtypes]
    for s in 1:nsites
        push!(sites[sitetypes[s]], s)
    end
    for b in 1:nbonds
        push!(bonds[bondtypes[b]], b)
    end
    Lattice(dim, sz, nsitetypes, nbondtypes,
            sites, bonds, nsites, nbonds, sitetypes, bondtypes,
            transvector, coords, bond_dirs, neighborsites, neighborbonds,
            source,target)
end
chain_lattice(param::Dict) = chain_lattice(param["L"])

@doc doc"""
    square_lattice(L::Integer, W::Integer=L)
    square_lattice(param::Dict)
    
Generates square lattice with size `L=param["L"]` $\times$ `W=param["W"]`.

`nsitetypes` is 2: all `1` (`2`) sites do not connect to another `1` (`2`) site.

`nbondtypes` is 2: `1` bonds are parallel to `x` axis and `2` are parallel to `y` axis.
"""
square_lattice(L::Integer) = square_lattice(L,L)
function square_lattice(L::Integer, W::Integer)
    s2xy(s::Integer) = mod(s-1,L),div(s-1,L)
    xy2s(x::Integer, y::Integer) = mod(y,W)*L + mod(x,L) + 1

    dim = 2
    sz = [L,W]
    nsitetypes = 2
    nbondtypes = 2
    nsites = L*W
    nbonds = 2*L*W
    sitetypes = zeros(Int, nsites)
    bondtypes = zeros(Int, nbonds)
    transvector = @compat Matrix(1.0I, dim, dim)
    coords = zeros(dim, nsites)
    bond_dirs = zeros(dim, nbonds)
    neighborsites = Vector{Int}[]
    neighborbonds = Vector{Int}[]
    source = zeros(Int,nbonds)
    target = zeros(Int,nbonds)
    sitetypes = zeros(Int,nsites)
    bondtypes = zeros(Int,nbonds)
    ns = zeros(Int,4)
    nb = zeros(Int,4)
    @inbounds for s in 1:nsites
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
        bond_dirs[:, 2s-1] = [1.0, 0.0]
        bond_dirs[:, 2s] = [0.0, 1.0]
        bondtypes[2s-1] = 1 
        bondtypes[2s] = 2 
    end
    sites = Vector{Int}[zeros(Int,0) for i in 1:nsitetypes]
    bonds = Vector{Int}[zeros(Int,0) for i in 1:nbondtypes]
    for s in 1:nsites
        push!(sites[sitetypes[s]], s)
    end
    for b in 1:nbonds
        push!(bonds[bondtypes[b]], b)
    end
    Lattice(dim, sz, nsitetypes, nbondtypes,
            sites, bonds, nsites, nbonds, sitetypes, bondtypes,
            transvector, coords, bond_dirs, neighborsites, neighborbonds,
            source,target)
end
function square_lattice(param::Dict)
    L = param["L"]
    W = get(param, "W", L)
    return square_lattice(L, W)
end

@doc doc"""
    triangular_lattice(L::Integer, W::Integer=L)
    triangular_lattice(param::Dict)
    
Generates triangular lattice with size `L=param["L"]` $\times$ `W=param["W"]`.

`nsitetypes` is 3: all `1`, `2`, and `3` sites do not connect to another `1`, `2`, and `3` site, respectively.

`nbondtypes` is 3: `1`, `2`, and `3` bonds make an angle of 0, 60, and 120 degree with `x` axis, respectively.
"""
triangular_lattice(L::Integer) = triangular_lattice(L,L)
function triangular_lattice(L::Integer, W::Integer)
    s2xy(s::Integer) = mod(s-1,L),div(s-1,L)
    xy2s(x::Integer, y::Integer) = mod(y,W)*L + mod(x,L) + 1
    ex = [1.0, 0.0]
    ey = [0.5, 0.5sqrt(3)]

    dim = 2
    sz = [L,W]
    nsitetypes = 3
    nbondtypes = 3
    nsites = L*W
    nbonds = 3nsites
    sitetypes = zeros(Int,nsites)
    bondtypes = zeros(Int,nbonds)
    transvector = [1.0 0.5; 0.0 0.5sqrt(3)]
    coords = zeros(dim, nsites)
    bond_dirs = zeros(dim, nbonds)
    neighborsites = Vector{Int}[]
    neighborbonds = Vector{Int}[]
    source = zeros(Int,nbonds)
    target = zeros(Int,nbonds)
    ns = zeros(Int,6)
    nb = zeros(Int,6)
    @inbounds for s in 1:nsites
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
        bond_dirs[:, 3s-2] .= [1.0,0.0]
        bond_dirs[:, 3s-1] .= [0.0,1.0]
        bond_dirs[:, 3s] .= [-1.0,1.0]
        bondtypes[3s-2] = 1
        bondtypes[3s-1] = 2
        bondtypes[3s] = 3
    end
    sites = Vector{Int}[zeros(Int,0) for i in 1:nsitetypes]
    bonds = Vector{Int}[zeros(Int,0) for i in 1:nbondtypes]
    for s in 1:nsites
        push!(sites[sitetypes[s]], s)
    end
    for b in 1:nbonds
        push!(bonds[bondtypes[b]], b)
    end
    Lattice(dim, sz, nsitetypes, nbondtypes,
            sites, bonds, nsites, nbonds, sitetypes, bondtypes,
            transvector, coords, bond_dirs, neighborsites, neighborbonds,
            source,target)
end
function triangular_lattice(param::Dict)
    L = param["L"]
    W = get(param, "W", L)
    return triangular_lattice(L, W)
end

@doc doc"""
    cubic_lattice(L::Integer, W::Integer=L, H::Integer=W)
    cubic_lattice(param::Dict)
    
Generates cubic lattice with size `L=param["L"]` $\times$ `W=param["W"]` $\times$ `H=param["H"]`.

`nsitetypes` is 2: all `1` (`2`) sites do not connect to another `1` (`2`) site.

`nbondtypes` is 2: `1`, `2`, and `3` bonds are parallel to `x`, `y`, and `z` axis, respectively.
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
    nsitetypes = 2
    nbondtypes = 3
    nsites = L*W*H
    nbonds = 3nsites
    sitetypes = zeros(Int,nsites)
    bondtypes = zeros(Int,nbonds)
    transvector = @compat Matrix(1.0I, dim, dim)
    coords = zeros(dim, nsites)
    bond_dirs = zeros(dim, nbonds)
    neighborsites = Vector{Int}[]
    neighborbonds = Vector{Int}[]
    source = zeros(Int,nbonds)
    target = zeros(Int,nbonds)
    ns = zeros(Int,6)
    nb = zeros(Int,6)
    @inbounds for s in 1:nsites
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
        bondtypes[3s-2] = 1
        bondtypes[3s-1] = 2
        bondtypes[3s] = 3
        bond_dirs[:, 3s-2] = [1.0, 0.0, 0.0]
        bond_dirs[:, 3s-1] = [0.0, 1.0, 0.0]
        bond_dirs[:, 3s-0] = [0.0, 0.0, 1.0]
    end
    sites = Vector{Int}[zeros(Int,0) for i in 1:nsitetypes]
    bonds = Vector{Int}[zeros(Int,0) for i in 1:nbondtypes]
    for s in 1:nsites
        push!(sites[sitetypes[s]], s)
    end
    for b in 1:nbonds
        push!(bonds[bondtypes[b]], b)
    end
    Lattice(dim, sz, nsitetypes, nbondtypes,
            sites, bonds, nsites, nbonds, sitetypes, bondtypes,
            transvector, coords, bond_dirs,
            neighborsites, neighborbonds,
            source,target)
end
function cubic_lattice(param::Dict)
    L = param["L"]
    W = get(param, "W", L)
    H = get(param, "H", W)
    return cubic_lattice(L, W, H)
end

@doc doc"""
    fully_connected_lattice(N::Integer)
    fully_connected_lattice(param::Dict)
    
Generates `N=param["N"]` site fully connected lattice

Both `nsitetypes` and `nbondtypes` are 1.
"""
function fully_connected_lattice(N::Integer)
    dim = 1
    nbonds = div(N*(N-1),2)
    nsitetypes = 1
    nbondtypes = 1
    sites = [collect(1:N)]
    bonds = [collect(1:nbonds)]
    sitetypes = ones(N)
    bondtypes = ones(nbonds)
    neighborsites = [vcat(1:s-1, s+1:N) for s in 1:N]
    neighborbonds = Vector{Int}[zeros(Int,0) for s in 1:N]
    source = zeros(Int,nbonds)
    target = zeros(Int,nbonds)
    ib = 0
    for i in 1:N
        for j in (i+1):N
            ib += 1
            source[ib] = i
            target[ib] = j
            push!(neighborbonds[i], ib)
            push!(neighborbonds[j], ib)
        end
    end

    lat = Lattice(dim, [N], nsitetypes, nbondtypes, sites, bonds, N, nbonds, sitetypes, bondtypes,
            zeros(1,1), zeros(1,N), zeros(1,nbonds), neighborsites, neighborbonds, source, target)
    return lat
end
fully_connected_lattice(param::Dict) = fully_connected_lattice(param["N"])
