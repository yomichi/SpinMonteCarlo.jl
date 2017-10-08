@testset "Lattice" begin

    param = Dict{String,Any}("L" => 3)
    @testset "$lname" for (lname, latt) in [("dimer", dimer_lattice),
                                            ("chain", chain_lattice),
                                            ("square", square_lattice),
                                            ("triangular", triangular_lattice),
                                            ("cubic", cubic_lattice),
                                           ]
        lat = latt(param)
        for s in 1:numsites(lat)
            cs = sitecoordinate(lat,s)
            for (ns, nb) in neighbors(lat, s)
                if source(lat,nb) == s
                    @test target(lat,nb) == ns
                    src = s
                    tgt = ns
                else
                    @test source(lat,nb) == ns
                    @test target(lat,nb) == s
                    src = ns
                    tgt = s
                end
                cs = lattice_sitecoordinate(lat,src)
                ct = lattice_sitecoordinate(lat,tgt)
                cb  = lattice_bonddirection(lat,nb)
                ct2 = (cs .+ cb .+ size(lat)) .% size(lat)
                @test ct == ct2
            end
        end
    end
end
