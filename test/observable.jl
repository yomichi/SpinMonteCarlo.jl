using Random

function test_single(scalartype, vectortype)
    nobs = 3
    ndata = 100
    obs_scalar = [scalartype() for i in 1:nobs]
    obs_vector = vectortype()
    rng = Random.MersenneTwister(SEED)
    X = randn(rng, nobs, ndata)
    for i in 1:ndata
        for j in 1:nobs
            push!(obs_scalar[j], X[j, i])
        end
        push!(obs_vector, X[:, i])
    end

    mean_scalar = [mean(obs_scalar[i]) for i in 1:nobs]
    mean_vector = mean(obs_vector)
    var_scalar = [var(obs_scalar[i]) for i in 1:nobs]
    var_vector = var(obs_vector)

    @test mean_scalar ≈ mean(X; dims=2)[:]
    @test var_scalar ≈ var(X; dims=2)[:]
    @test mean_vector ≈ mean(X; dims=2)[:]
    @test var_vector ≈ var(X; dims=2)[:]
end

function test_binning(scalartype, vectortype)
    nobs = 3
    ndata = 100
    obs_scalar = [scalartype() for i in 1:nobs]
    obs_vector = vectortype()
    rng = Random.MersenneTwister(SEED)
    X = randn(rng, nobs, ndata)
    for i in 1:ndata
        for j in 1:nobs
            push!(obs_scalar[j], X[j, i])
        end
        push!(obs_vector, X[:, i])
    end

    binning_scalar = [binning(obs_scalar[i]) for i in 1:nobs]
    binning_vector = binning(obs_vector)

    mean_scalar = [mean(binning_scalar[i]) for i in 1:nobs]
    mean_vector = mean(binning_vector)
    var_scalar = [var(binning_scalar[i]) for i in 1:nobs]
    var_vector = var(binning_vector)

    @test mean_scalar ≈ mean_vector
    @test var_scalar ≈ var_vector
end

function test_jackknife(scalartype, vectortype)
    nobs = 3
    ndata = 100
    obs_scalar_x = [scalartype() for i in 1:nobs]
    obs_vector_x = vectortype()
    obs_scalar_y = [scalartype() for i in 1:nobs]
    obs_vector_y = vectortype()
    rng = Random.MersenneTwister(SEED)
    X = randn(rng, nobs, ndata)
    Y = randn(rng, nobs, ndata)
    for i in 1:ndata
        for j in 1:nobs
            push!(obs_scalar_x[j], X[j, i])
            push!(obs_scalar_y[j], Y[j, i])
        end
        push!(obs_vector_x, X[:, i])
        push!(obs_vector_y, Y[:, i])
    end
    jk_scalar_x = [jackknife(obs_scalar_x[i]) for i in 1:nobs]
    jk_scalar_y = [jackknife(obs_scalar_y[i]) for i in 1:nobs]

    sxpcy_scalar = [sin(jk_scalar_x[i]) + cos(jk_scalar_y[i]) for i in 1:nobs]
    means_scalar = [mean(sxpcy_scalar[i]) for i in 1:nobs]
    vars_scalar = [var(sxpcy_scalar[i]) for i in 1:nobs]

    jk_vector_x = jackknife(obs_vector_x)
    jk_vector_y = jackknife(obs_vector_y)

    sxpcy_vector = sin(jk_vector_x) + cos(jk_vector_y)
    means_vector = mean(sxpcy_vector)
    vars_vector = var(sxpcy_vector)

    @test means_scalar ≈ means_vector
    # @test vars_scalar ≈ vars_vector

    jk_scalar_x[1] + jk_vector_x
    return jk_scalar_x[1] * jk_vector_x
end

@testset "Simple" begin
    @testset "Simple" begin
        test_single(SimpleObservable, SimpleVectorObservable)
    end
    @testset "Tiny" begin
        test_single(TinyObservable, TinyVectorObservable)
    end
end
@testset "Binning" begin
    @testset "Simple" begin
        test_binning(SimpleObservable, SimpleVectorObservable)
    end
end
@testset "Jackknife" begin
    @testset "Simple" begin
        test_jackknife(SimpleObservable, SimpleVectorObservable)
    end
end

@testset "Observable bug fixes" begin
    @testset "binning validates inputs" begin
        scalar = SimpleObservable()
        vector = SimpleVectorObservable()
        @test_throws ArgumentError binning(scalar)
        @test_throws ArgumentError binning(vector)

        push!(scalar, 1.0)
        push!(scalar, 2.0)
        push!(vector, [1.0, 2.0])
        push!(vector, [3.0, 4.0])
        @test_throws ArgumentError binning(scalar; numbins=3)
        @test_throws ArgumentError binning(vector; numbins=3)

        param = Parameter("Model" => Ising,
                          "Lattice" => "fully connected graph",
                          "N" => 4,
                          "T" => 1.0,
                          "Update Method" => local_update!,
                          "MCS" => 16,
                          "Thermalization" => 0,
                          "Seed" => SEED)
        @test runMC(param) isa AbstractDict
    end

    @testset "SimpleVectorObservable merge" begin
        lhs = SimpleVectorObservable()
        rhs = SimpleVectorObservable()
        push!(lhs, [1.0, 2.0])
        push!(rhs, [3.0, 5.0])

        merge!(lhs, rhs)
        @test count(lhs) == 2
        @test mean(lhs) ≈ [2.0, 3.5]
        @test var(lhs) ≈ [2.0, 4.5]
    end

    @testset "JackknifeVector var initializes accumulator" begin
        jk = JackknifeVector([[1.0, 2.0], [3.0, 5.0], [5.0, 8.0]])
        @test var(jk) ≈ [16.0 / 3.0, 12.0]
    end

    @testset "JackknifeVector manual broadcast methods" begin
        jk = JackknifeVector([[1.0, 2.0], [3.0, 4.0]])
        @test mean(broadcast(*, jk, 2.0)) ≈ [4.0, 6.0]
    end

    @testset "jackknife vector observable set" begin
        oset = SimpleVectorObservableSet()
        oset["x"] = SimpleVectorObservable()
        push!(oset["x"], [1.0, 2.0])
        push!(oset["x"], [3.0, 4.0])

        jkset = jackknife(oset)
        @test jkset isa JackknifeVectorSet
        @test jkset["x"] isa JackknifeVector
        @test mean(jkset["x"]) ≈ [2.0, 3.0]
    end

    @testset "BinningVectorObservable rebinning" begin
        nobs = 3
        ndata = 32
        minbinnum = 4
        rng = Random.MersenneTwister(SEED)
        X = randn(rng, nobs, ndata)

        vector = BinningVectorObservable(minbinnum)
        scalars = [BinningObservable(minbinnum) for j in 1:nobs]
        for i in 1:ndata
            push!(vector, X[:, i])
            for j in 1:nobs
                push!(scalars[j], X[j, i])
            end
        end

        @test SpinMonteCarlo.maxlevel(vector) == SpinMonteCarlo.maxlevel(scalars[1])
        @test count(vector) == ndata
        @test mean(vector) ≈ [mean(s) for s in scalars]
        @test var(vector) ≈ [var(s) for s in scalars]
        @test stderror(vector) ≈ [stderror(s) for s in scalars]
        @test tau(vector) ≈ [tau(s) for s in scalars]
    end

    @testset "BinningVectorObservable extrapolation" begin
        nobs = 3
        ndata = 1024
        minbinnum = 4
        rng = Random.MersenneTwister(SEED)
        X = randn(rng, nobs, ndata)

        vector = BinningVectorObservable(minbinnum)
        scalars = [BinningObservable(minbinnum) for j in 1:nobs]
        for i in 1:ndata
            push!(vector, X[:, i])
            for j in 1:nobs
                push!(scalars[j], X[j, i])
            end
        end

        for extrapolate in (extrapolate_stderror, extrapolate_tau)
            values, errors = extrapolate(vector)
            reference = [extrapolate(s) for s in scalars]
            @test values ≈ [r[1] for r in reference]
            @test errors ≈ [r[2] for r in reference]
        end
    end

    @testset "zeros of observables" begin
        ## the dimensions can be given either as a vararg or as a tuple
        for dims in ((2, 3), ((2, 3),))
            obs = zeros(SimpleObservable, dims...)
            @test size(obs) == (2, 3)
            @test eltype(obs) === SimpleObservable

            ## each element has to be an independent object
            push!(obs[1, 1], 1.0)
            @test count(obs[1, 1]) == 1
            @test count(obs[2, 3]) == 0
        end
        @test size(zeros(SimpleVectorObservable, 2)) == (2,)
        @test size(zeros(TinyVectorObservable, (2,))) == (2,)
    end

    @testset "observable set helpers" begin
        oset = SimpleObservableSet()
        makeMCObservable!(oset, "x")
        push!(oset["x"], 1.0)

        @test_logs (:warn, r"already exists") makeMCObservable!(oset, "x")
        reset!(oset)
        @test count(oset["x"]) == 0
    end
end
