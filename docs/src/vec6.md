```@meta
CurrentModule = ExpMap
```

# VEC6 — Immutable 6D Vector

`VEC6` is an **immutable** struct representing a 6-component vector with
`Float64` entries indexed `x1` … `x6`.

In rigid-body mechanics a `VEC6` typically encodes a **spatial** (6D) quantity:

| Usage     | Layout                        |
|-----------|-------------------------------|
| Wrench    | `[force(3); torque(3)]`       |
| Twist     | `[linear vel.(3); angular vel.(3)]` |
| State     | `[position(3); rotation vec.(3)]`   |

!!! note "Immutability"
    `setindex` (without `!`) returns a **new** `VEC6`; the original is unchanged.
    Scalar multiplication accepts any `Real`, so integer literals like `2 * v` work.

## VEC6 construction

| Syntax                              | Description                                     |
|-------------------------------------|-------------------------------------------------|
| `VEC6(x1,x2,x3,x4,x5,x6)`          | From six `Float64` scalars                      |
| `VEC6(t::NTuple{6,Float64})`        | From a 6-tuple                                  |
| `VEC6(v::Vector{Float64})`          | From a 6-element Julia vector (asserts length)  |
| `VEC6(a::VEC3, b::VEC3)`            | Concatenate two `VEC3`: `[a; b]`                |
| `VEC6(a::VEC3, b::RV3)`             | Concatenate a `VEC3` and an `RV3`: `[a; b]`     |
| `zero(VEC6)`                        | Zero vector `(0,0,0,0,0,0)`                     |
| `VEC6_unit(k, a)`                   | Component `k` set to `a`, rest zero             |

## VEC6 arithmetic operators

| Expression         | Result  | Description                       |
|--------------------|---------|-----------------------------------|
| `a + b`            | `VEC6`  | Component-wise addition           |
| `a - b`            | `VEC6`  | Component-wise subtraction        |
| `-v`               | `VEC6`  | Unary negation                    |
| `λ * v`, `v * λ`  | `VEC6`  | Scalar multiplication (`λ::Real`) |
| `v / λ`            | `VEC6`  | Scalar division                   |

## VEC6 geometric operations

| Function    | Returns   | Description         |
|-------------|-----------|---------------------|
| `norm2(v)`  | `Float64` | Euclidean norm ‖v‖  |

## VEC6 accessors and indexing

| Syntax                   | Description                                          |
|--------------------------|------------------------------------------------------|
| `v.x1` … `v.x6`          | Direct field access                                  |
| `v.v`                    | Convert to `Vector{Float64}` `[x1,…,x6]`             |
| `v[i]`                   | 1-based integer indexing                             |
| `setindex(v, val, i)`    | New `VEC6` with component `i` set to `val`           |
| `copy(v)`                | Deep copy                                            |

```@docs
VEC6
VEC6_unit
norm2(v::VEC6)
```
