abstract Model

type Ising <: Model
    lat :: Lattice
    spins :: Vector{Int}

    function Ising(lat::Lattice)
        model = new()
        model.lat = lat
        model.spins = rand([1,-1], numsites(lat))
        return model
    end
end

type Potts <: Model
    lat :: Lattice
    Q :: Int
    spins :: Vector{Int}

    function Potts(lat::Lattice, Q::Integer)
        spins = rand(1:Q, numsites(lat))
        return new(lat, Q, spins)
    end
end

type Clock <: Model
    lat :: Lattice
    Q :: Int
    spins :: Vector{Int}
    cosines :: Vector{Float64}
    sines :: Vector{Float64}

    function Clock(lat::Lattice, Q::Integer)
        spins = rand(1:Q, numsites(lat))
        cosines = [cospi(2s/Q) for s in 1:Q]
        sines = [sinpi(2s/Q) for s in 1:Q]
        return new(lat, Q, spins, cosines, sines)
    end
end

type XY <: Model
    lat :: Lattice
    spins :: Vector{Float64}

    function XY(lat::Lattice)
        model = new()
        model.lat = lat
        model.spins = rand(numsites(lat))
        return model
    end
end

