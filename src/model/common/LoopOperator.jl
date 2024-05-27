@doc doc"""
Enumtype including `LET_*`
"""
@enum(LoopElementType,
      LET_Cut,    # [1 1; 1 1]
      LET_FMLink, # [1 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 1]
      LET_AFLink, # [0 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]
      LET_Vertex, # [0 0 0 0; 0 1 1 0; 0 1 1 0; 0 0 0 0]
      LET_Cross,)

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
    let_type::LoopElementType
    isdiagonal::Bool
    time::Float64
    space::Int
    subspace::Tuple{Int,Int}
    bottom_id::Int
    top_id::Int
end
function LocalLoopOperator(let_type::LoopElementType, time::Real, space::Int,
                           subspace::Tuple{Int,Int})
    return LocalLoopOperator(let_type, true, time, space, subspace, 0, 0)
end
