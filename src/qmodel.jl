@doc doc"""
Enumtype including `LET_*`
"""
@enum(LoopElementType,
      LET_Cut,    # [1 1; 1 1]
      LET_FMLink, # [1 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 1]
      LET_AFLink, # [0 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]
      LET_Vertex, # [0 0 0 0; 0 1 1 0; 0 1 1 0; 0 0 0 0]
      LET_Cross,  # [1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1]
     )

@doc """
Loop element depicted as

```
|
o
|
```

or matrix 

```
1 1 |+>
1 1 |->
```
""" LET_Cut

@doc """
Loop element depicted as

```
|  |
|--|
|  |
```

or matrix 

```
1 0 0 0 |++>
0 0 0 0 |+->
0 0 0 0 |-+>
0 0 0 1 |-->
```
""" LET_FMLink

@doc """
Loop element depicted as

```
|  |
|~~|
|  |
```

or matrix 

```
0 0 0 0 |++>
0 1 0 0 |+->
0 0 1 0 |-+>
0 0 0 0 |-->
```
""" LET_AFLink

@doc """
Loop element depicted as

```
|  |
~~~~
~~~~
|  |
```

or matrix 

```
0 0 0 0 |++>
0 1 1 0 |+->
0 1 1 0 |-+>
0 0 0 0 |-->
```
""" LET_Vertex

@doc raw"""
Loop element depicted as

```
|  |
 \/
 /\
|  |
```

or matrix 

```
1 0 0 0 |++>
0 0 1 0 |+->
0 1 0 0 |-+>
0 0 0 1 |-->
```
""" LET_Cross

@doc doc"""
(Imaginary-temporary and spatial) local operator as a perturbation with assigned loop element.
# Fields
- `let_type` : assigned loop element
- `isdiagonal` : operator is diagonal or not
    - in other words, two states connecting this perturbation are equivalent to
      each other or not.
- `time` : imaginary time ($\tau/\beta \in [0,1)$) which this perturbation acts on.
- `space` : spin or bond which this perturbation acts on.
            denotes `space` spin if space <= nspins
                or  `space - nspins` bond otherwise.
- `subspace` : subspin(s) index
- `bottom_id` :: index of node of union find assigned to a loop
- `top_id` :: index of node of union find assigned to the other loop

"""
mutable struct LocalLoopOperator
    let_type :: LoopElementType
    isdiagonal :: Bool
    time :: Float64
    space :: Int
    subspace :: Tuple{Int,Int}
    bottom_id :: Int
    top_id :: Int
end
function LocalLoopOperator(let_type::LoopElementType, time::Real, space::Int, subspace::Tuple{Int,Int})
    LocalLoopOperator(let_type, true, time, space, subspace, 0,0)
end

abstract type QuantumLocalZ2Model <: Model end

@doc doc"""
Spin-$S$ XXZ model represented as the following Hamiltonian,
\begin{equation}
\mathcal{H} = \sum_{i,j} \left[ J_{ij}^z S_i^z S_j^z 
            + \frac{J_{ij}^{xy}}{2} (S_i^+ S_j^- + S_i^-S_j^+) \right]
            - \sum_i \Gamma_i S_i^x,
\end{equation}
where $S^x, S^y, S^z$ are $x, y$ and $z$ component of spin operator with length $S$,
and $S^\pm \equiv S^x \pm iS^y$ are ladder operator.
A state is represented by a product state (spins at $\tau=0$) of local $S^z$ diagonal basis and an operator string (perturbations).
A local spin with length $S$ is represented by a symmetrical summation of $2S$ sub spins with length $1/2$.
"""
mutable struct QuantumXXZ <: QuantumLocalZ2Model
    lat :: Lattice
    S2 :: Int
    spins :: Vector{Int}
    ops :: Vector{LocalLoopOperator}
    rng :: Random.MersenneTwister
    parameter :: Parameter

    function QuantumXXZ(lat::Lattice, S::Real)
        if round(2S) != 2S
            error("`S` should be integer or half-integer")
        end
        S2 = round(Int,2S)
        model = new()
        model.rng = seed!(Random.MersenneTwister(0))
        model.lat = lat
        model.S2 = S2
        model.spins = rand(model.rng,[1,-1], numsites(lat)*S2)
        model.ops = LocalLoopOperator[]
        model.parameter = Parameter()
        return model
    end
    function QuantumXXZ(lat::Lattice, S::Real, seed)
        if round(2S) != 2S
            error("`S` should be integer or half-integer")
        end
        S2 = round(Int,2S)
        model = new()
        model.rng = seed!(Random.MersenneTwister(0), seed)
        model.lat = lat
        model.S2 = S2
        model.spins = rand(model.rng,[1,-1], numsites(lat)*S2)
        model.ops = LocalLoopOperator[]
        model.parameter = Parameter()
        return model
    end
end
@doc doc"""
    QuantumXXZ(param)

Generates `QuantumXXZ` using `param["Lattice"]`, `param["S"]` and `param["Seed"]` (if defined).
Each subspin will be initialized independently and randomly.
"""
function QuantumXXZ(param::Parameter, lat::Lattice)
    S = param["S"]
    if "Seed" in keys(param)
        model = QuantumXXZ(lat,S,param["Seed"])
    else
        model = QuantumXXZ(lat,S)
    end
    merge!(model.parameter, param)
    return model
end
QuantumXXZ(param::Parameter) = QuantumXXZ(param, generatelattice(param))

@inline site2subspin(site::Integer, ss::Integer, S2::Integer) = (site-1)*S2+ss
@inline function subspin2site(subspin::Integer, S2::Integer)
    return ceil(Int, subspin/S2), mod1(subspin,S2)
end
