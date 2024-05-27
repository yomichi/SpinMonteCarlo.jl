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
