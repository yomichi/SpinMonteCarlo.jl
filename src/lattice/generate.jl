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

function generatelattice(param
             ; 
             bravaisdict::Dict{String,LatticeParameter}=stdbravais,
             unitcelldict::Dict{String, LatticeParameter}=stdunitcells,
             latticedict::Dict{String, LatticeParameter}=stdlattices
            )
    lat = latticedict[param["Lattice"]]
    bra = bravaisdict[lat.bravais]
    ucell = unitcelldict[lat.unitcell]
    @assert lat.dimension == bra.dimension == ucell.dimension
    D = lat.dimension
    bparams = Dict(name => get(param,String(name),default) for (name, default) in bra.parameters)

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

    sites = ucell.sites
    numsites_in_cell = length(sites)
    numsites = numcell * numsites_in_cell

    let x = Int[]
        for site in sites
            push!(x, site.id)
        end
        sort!(x)
        @assert x == collect(1:numsites_in_cell)
    end

    bonds = ucell.bonds

    g = MetaGraph(numsites)

    numbonds = 0
    numsitetypes = 0
    numbondtypes = 0
    nss = Int[]
    nbs = Int[]

    for icell in 0:(numcell-1)
        cellcoord = index2coord(icell, L)
        for site in sites
            id = numsites_in_cell * icell + site.id
            set_prop!(g, id, :sitetype, site.sitetype)
            while site.sitetype > numsitetypes
                numsitetypes += 1
                push!(nss, 0)
            end
            nss[site.sitetype] += 1
            coord = latvec * (cellcoord .+ site.coord)
            set_prop!(g, id, :coord, coord)
            set_prop!(g, id, :localsite, site.id)
            set_prop!(g, id, :cellcoord, cellcoord)
            for bond in bonds
                source = numsites_in_cell * coord2index(cellcoord .+ bond.source.offset, L) + bond.source.id
                target = numsites_in_cell * coord2index(cellcoord .+ bond.target.offset, L) + bond.target.id
                add_edge!(g, source, target)
                numbonds += 1
                set_prop!(g, source, target, :bondtype, bond.bondtype)
                set_prop!(g, source, target, :source, source)
                set_prop!(g, source, target, :target, target)
                while bond.bondtype > numbondtypes
                    numbondtypes += 1
                    push!(nbs, 0)
                end
                nbs[bond.bondtype] += 1
            end
        end
    end

    set_prop!(g, :dim, D)
    set_prop!(g, :numsites, numsites)
    set_prop!(g, :numbonds, numbonds)
    set_prop!(g, :size, L)
    set_prop!(g, :numsitetypes, numsitetypes)
    set_prop!(g, :numbondtypes, numbondtypes)
    set_prop!(g, :nss, nss)
    set_prop!(g, :nbs, nbs)

    for bond in edges(g)
        src = sitecoordinate(g, source(g, bond))
        tgt = sitecoordinate(g, target(g, bond))
        dir = rem.(tgt .- src, latvec*L, RoundNearest)
        set_prop!(g, bond, :bonddirection, dir)
        set_prop!(g, bond, :weight, norm(dir))
    end
    return g
end
