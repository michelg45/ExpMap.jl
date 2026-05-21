```@meta
CurrentModule = ExpMap
```

# VEC3 — Immutable 3D Vector

`VEC3` is a lightweight, **immutable** struct representing a 3D vector with
`Float64` components `x`, `y`, `z`.  It is the fundamental spatial primitive
used throughout ExpMap.

## VEC3 construction

| Syntax                          | Description                              |
|---------------------------------|------------------------------------------|
| `VEC3(x, y, z)`                 | From three `Float64` scalars             |
| `VEC3(v::Vector{Float64})`      | From a 3-element Julia vector            |
| `VEC3()`                        | Zero vector `(0.0, 0.0, 0.0)`            |
| `VEC3(r::RV3)`                  | Reinterpret an `RV3` as a plain `VEC3`   |

## VEC3 arithmetic operators

All operators are **component-wise** unless stated otherwise.

| Expression      | Result  | Description                          |
|-----------------|---------|--------------------------------------|
| `v1 + v2`       | `VEC3`  | Addition                             |
| `v1 - v2`       | `VEC3`  | Subtraction                          |
| `-v`            | `VEC3`  | Unary negation                       |
| `v1 * v2`       | `VEC3`  | Hadamard (component-wise) product    |
| `v * s`, `s * v`| `VEC3`  | Scalar multiplication                |
| `v1 / v2`       | `VEC3`  | Component-wise division              |
| `v / s`         | `VEC3`  | Scalar division                      |
| `s / v`         | `VEC3`  | Scalar divided by each component     |

!!! note "Hadamard vs. geometric products"
    `v1 * v2` is the **component-wise** (Hadamard) product — not the dot or
    cross product.  Use `dotp(v1, v2)` and `crossp(v1, v2)` for those.

## VEC3 geometric operations

| Function            | Returns   | Description                   |
|---------------------|-----------|-------------------------------|
| `dotp(v1, v2)`      | `Float64` | Dot product  v1 · v2          |
| `crossp(v1, v2)`    | `VEC3`    | Cross product  v1 × v2        |
| `norm2(v)`          | `Float64` | Euclidean norm  ‖v‖           |

## VEC3 accessors and indexing

| Syntax       | Description                                         |
|--------------|-----------------------------------------------------|
| `v.x`, `v.y`, `v.z` | Direct field access                       |
| `v.v`        | Convert to `Vector{Float64}` `[x, y, z]`           |
| `v[1]`…`v[3]`| 1-based integer indexing (1→x, 2→y, 3→z)          |
| `copy(v)`    | Deep copy                                           |

```@docs
VEC3
dotp(v1::VEC3, v2::VEC3)
crossp(v1::VEC3, v2::VEC3)
norm2(v::VEC3)
```
