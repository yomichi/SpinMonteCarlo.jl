using SpinMonteCarlo
using Makie

latplot(lat::Lattice) = latplot(lat, Val{dim(lat)})

function latplot(lat::Lattice, ::Type{Val{1}})
    points = Point2f0[]
    minx = Inf
    maxx = -Inf
    for site in sites(lat)
        c = sitecoordinate(site)[1]
        push!(points, Point2(c,0))
        minx = min(minx, c)
        maxx = max(maxx, c)
    end
    limits = FRect(minx-1.0, -1.0, maxx-minx+2.0, 2.0)
    segments = Pair{Point2f0, Point2f0}[]
    for bond in bonds(lat)
        src = sitecoordinate(lat, source(bond))[1]
        tgt = bonddirection(bond)[1] + src
        push!(segments, Point2f0(src,0)=>Point2f0(tgt,0))
    end
    scene = scatter(points, limits=limits)
    linesegments!(scene, segments, limits=limits)
    axis = scene[Axis]
    axis[:names, :axisnames] = ("", "")
    axis[:frame, :linewidth] = 0.0
    axis[:grid, :linewidth] = (0.0, 0.0)
    scene
end

function latplot(lat::Lattice, ::Type{Val{2}})
    points = Point2f0[]
    minx = miny = Inf
    maxx = maxy = -Inf
    for site in sites(lat)
        c = sitecoordinate(site)
        push!(points, Point2f0(c))
        minx = min(minx, c[1])
        maxx = max(maxx, c[1])
        miny = min(miny, c[2])
        maxy = max(maxy, c[2])
    end
    limits = FRect(minx-1.0, miny-1.0, maxx-minx+2.0, maxy-miny+2.0)
    segments = Pair{Point2f0, Point2f0}[]
    for bond in bonds(lat)
        src = sitecoordinate(lat, source(bond))
        tgt = bonddirection(bond) .+ src
        push!(segments, Point2f0(src)=>Point2f0(tgt))
    end
    scene = scatter(points, limits=limits)
    linesegments!(scene, segments, limits=limits)
    axis = scene[Axis]
    axis[:names, :axisnames] = ("", "")
    axis[:frame, :linewidth] = 0.0
    axis[:grid, :linewidth] = (0.0, 0.0)
    scene
end

function latplot(lat::Lattice, ::Type{Val{3}})
    points = Point3f0[]
    minx = miny = minz = Inf
    maxx = maxy = maxz = -Inf
    for site in sites(lat)
        c = sitecoordinate(site)
        push!(points, Point3f0(c))
        minx = min(minx, c[1])
        maxx = max(maxx, c[1])
        miny = min(miny, c[2])
        maxy = max(maxy, c[2])
        minz = min(minz, c[3])
        maxz = max(maxz, c[3])
    end
    segments = Pair{Point3f0, Point3f0}[]
    for bond in bonds(lat)
        src = sitecoordinate(lat, source(bond))
        tgt = bonddirection(bond) .+ src
        push!(segments, Point3f0(src)=>Point3f0(tgt))
    end
    scene = scatter(points)
    linesegments!(scene, segments)
    axis = scene[Axis]
    axis[:names, :axisnames] = ("", "", "")
    axis[:frame, :linewidth] = 0.0
    scene
end
