export generatelattice

function interpolate!(ex,env)
    if isa(ex,Symbol)
        if haskey(env, ex)
            return env[ex]
        end
    elseif isa(ex,Expr)
        for (i,x) in enumerate(ex.args)
            if isa(x,Symbol)
                if haskey(env, x)
                    ex.args[i] = env[x]
                end
            elseif isa(x,Expr)
                ex.args[i] = interpolate!(ex.args[i], env)
            end
        end
    end
    return ex
end

function index2coord(index::Integer, L::AbstractArray)
    D = length(L)
    coord = zeros(Int,D)
    for d in 1:D
        coord[d] = index%L[d]
        index ÷= L[d]
    end
    return coord
end

function coord2index(coord::AbstractArray, L::AbstractArray)
    D = length(L)
    index = 0
    offset = 0
    c = mod.(coord, L)
    for d in D:-1:1
        index *= L[d]
        index += c[d]
    end
    return index
end

@doc """
    generatelattice(param::Parameter)

generates `Lattice` from `Parameter`.
"""
function generatelattice(param)
    latticedict = get(param, "LatticeDict", stdlattices)
    lat = latticedict[param["Lattice"]]
    return generator(lat)(param)
end

function generatelattice_std(param)
    bravaisdict = get(param, "BravaisDict", stdbravais)
    unitcelldict = get(param, "UnitcellDict", stdunitcells)
    latticedict = get(param, "LatticeDict", stdlattices)
    lat = latticedict[param["Lattice"]]
    bra = bravaisdict[lat.bravais]
    ucell = unitcelldict[lat.unitcell]
    @assert lat.dimension == bra.dimension == ucell.dimension
    D = lat.dimension
    bparams = Dict(name => get(param,String(name),default) for (name, default) in bra.parameters)
    for (name, default) in lat.parameters
        bparams[name] = get(param,String(name),default)
    end

    bc = get(param, "Periodic Boudary Condition", lat.periodic)

    if isa(param["L"], Vector)
        L = param["L"]
    else
        L = [param["L"]]
    end
    while length(L) < D
        push!(L,L[end])
    end

    numcell = prod(L)

    latvec = interpolate!(bra.basis, bparams) |> eval

    usites = ucell.sites
    numsites_in_cell = length(usites)
    numsites = numcell * numsites_in_cell

    let x = Int[]
        for site in usites
            push!(x, site.id)
        end
        sort!(x)
        @assert x == collect(1:numsites_in_cell)
    end

    ubonds = ucell.bonds

    numbonds = 0
    numsitetypes = 0
    numbondtypes = 0
    nss = Int[]
    nbs = Int[]

    sites = Site[]
    bonds = Bond[]
    ib = 0

    use_index_as_sitetype = get(param, "Use Indecies as Site Types", false)
    use_index_as_bondtype = get(param, "Use Indecies as Bond Types", false)

    for icell in 0:(numcell-1)
        cellcoord = index2coord(icell, L)
        for site in usites
            id = numsites_in_cell * icell + site.id
            coord = latvec * (cellcoord .+ site.coord)
            s = Site(id, ifelse(use_index_as_sitetype, id, site.sitetype), Int[], Int[], coord, site.id, cellcoord)
            push!(sites,s)
        end
        for bond in ubonds
            for d in 1:D
                if !(bc[d] || bond.source.offset[d] == bond.source.offset[d])
                    if !(0 <= cellcoord[d] + (bond.target.offset[d] - bond.source.offset[d]) < L[d] )
                        continue
                    end
                end
            end
            source = numsites_in_cell * coord2index(cellcoord .+ bond.source.offset, L) + bond.source.id
            target = numsites_in_cell * coord2index(cellcoord .+ bond.target.offset, L) + bond.target.id
            dir = latvec * ((bond.target.offset .+ usites[bond.target.id].coord ) 
                            .- (bond.source.offset .+ usites[bond.source.id].coord ))
            ib += 1
            b = Bond(ib, ifelse(use_index_as_bondtype, ib, bond.bondtype), source, target, dir) 
            push!(bonds, b)
        end
    end

    return Lattice(latvec, L, sites, bonds)
end

function generate_fully_connected_graph(param)
    N = param["N"]
    nb = div(N*(N-1),2)
    sites = Site[]
    bonds = Bond[]
    ib = 0
    for s in 1:N
        push!(sites, Site(s, 1, Int[], Int[], [0.0], 1, [0]))
        for t in (s+1):N
            ib += 1
            push!(bonds, Bond(ib, 1, s, t, [0.0]))
        end
    end
    latvec = zeros(1,1)
    return Lattice(latvec, [N], sites, bonds)
end
