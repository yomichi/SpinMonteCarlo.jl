@testset "removed snapshot API" begin
    function removed_error(f, name)
        try
            f()
        catch err
            return err isa ErrorException &&
                   occursin("$(name) was removed in v1.3", err.msg) &&
                   occursin("A visualization extension is planned", err.msg)
        end
        return false
    end

    @test_throws ErrorException gen_snapshot!()
    @test_throws ErrorException gen_snapshot!(nothing, 1.0; MCS=1)
    @test removed_error(() -> gen_snapshot!(), "gen_snapshot!")

    @test_throws ErrorException gensave_snapshot!()
    @test_throws ErrorException gensave_snapshot!(IOBuffer(), nothing, 1.0)
    @test removed_error(() -> gensave_snapshot!(), "gensave_snapshot!")

    @test_throws ErrorException load_snapshot()
    @test_throws ErrorException load_snapshot(IOBuffer("1 2 3"))
    @test removed_error(() -> load_snapshot(), "load_snapshot")
end
