@testset "Lattice" begin

    param = Dict{String,Any}("L" => 3)
    @testset "$lname" for lname in ["dimer",
                                    "chain lattice",
                                    "square lattice",
                                    "triangular",
                                    "cubic lattice",
                                   ]
        lat = generatelattice(param)
        for s in 1:numsites(lat)
            cs = sitecoordinate(lat,s)
        end
    end
end
