using JSON

@testset "$modelname" for modelname in ("QuantumXXZ",)
    files = filter(s->endswith(s,".json"), readdir(joinpath("ref", modelname)))
    @testset "$filename" for filename in files
        diagres = JSON.parsefile(joinpath("ref",modelname,filename))
        param = diagres["Parameter"]
        ref = diagres["Result"]
        param["Model"] = QuantumXXZ
        param["Lattice"] = "chain lattice"
        param["Update Method"] = loop_update!

        mcres = runMC(param)
        @testset "$name" for name in keys(ref)
            p = p_value(mcres[name],ref[name])
            @test p > alpha
            if p <= alpha
                @show filename, name, mcres[name], ref[name]
            end
        end
    end
end

