# Lattice

`generatelattice(param::Parameter)` generates a `Lattice` from the parameter.

- "Lattice"
    - The name of the lattice
- "L"
    - Length along each dimension

If you want to use a lattice/bravairs/unitcell you defined, 
the following optional parameters must be set:

- "LatticeDict"
- "BravaisDict"
- "UnitcellDict"

## Pre-defined lattices

- "chain lattice"
- "bond-alternating chain lattice"
- "square lattice"
- "triangular lattice"
- "honeycomb lattice"
- "ladder"
- "cubic lattice"
- "fully connected lattice"

## Define your lattice
A lattice can be defined by combining a Bravais lattice (lattice basic vector) and a unit cell (sublattice structure).
A lattice, a Bravais lattice, and a unit cell are represented by an instance of `P = Dict{Symbol, Any}`.

### Lattice

- `:name => String`
- `:dimenstion => Integer`
- `:bravais => String`
- `:unitcell => String`
- `:parameters => Vector{Tuple{Symbol, Any}}`
    - To overwrite the default value of parameters of a Bravais lattice
    - An element of a value is a tuple of `parameter_name :: Symbol` and `new_default_value :: Any`.
    - ex. `[(:a, sqrt(3.0))]`
- `:periodic => Vector{Bool}`
    - To specify the boundary condition (BC) for each dimension, periodic BC if `true` and open BC if `false`.

### Bravais lattice

- `:name => String`
- `:dimension => Integer`
- `:parameters => Vector{Tuple{Symbol, Any}}`
    - parameter of the lattice vector
    - An element of a value is a tuple of `parameter_name :: Symbol` and `default_value :: Any`.
    - ex. `[(:a, 1.0)]`
- `:basis => Expr`
    - parameterized lattice basic vector
    - Type of a value is `Expr` of an `dim×dim Matrix`
    - ex. `:( [a 0.0; 0.0 a] )`

### Unit cell

- `:name => String`
- `:dimension => Integer`
- `:sites => Vector{P}`
- `:bonds => Vector{P}`

A site is an instance of `Dict{Symbol, Any}` having the following items:

- `:id => Integer`
    - The site index in a unitcell (1-origin)
- `:sitetype => Integer`
    - The site type (1-origin)
- `:coord => Vector{Float64}`
    - The fractional coordinate of the site

A bond is an instance of `Dict{Symbol, Any}` having the following items:

- `:bondtype => Integer`
    - The bond type (1-origin)
- `:source => P`
    - One end (site) of a bond
    - The value is a `P` having the following items:
        - `:id => Integer` : the site index in a unitcell
        - `:offset => Vector{Integer}` : the coordinate of the unitcell
- `:target => P`
    - Another end (site) of a bond
