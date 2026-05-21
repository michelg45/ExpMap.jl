```@meta
CurrentModule = ExpMap
```

# RV3 — Cartesian Rotation Vector

`RV3` is an **immutable** struct encoding a 3D rotation as a
**Cartesian Rotation Vector** (CRV), also known as a rotation pseudo-vector
or Euler vector.

A CRV encodes a rotation as  **ψ = θ · n̂**, where:
- **θ = ‖ψ‖** is the rotation angle (radians)
- **n̂ = ψ / ‖ψ‖** is the unit rotation axis

The zero vector `RV3()` represents the identity (no rotation).

!!! note "Immutability"
    `RV3` is immutable.  Operations that would modify a component return
    a **new** `RV3`; see `setindex`.

## RV3 construction

| Syntax                         | Description                                    |
|--------------------------------|------------------------------------------------|
| `RV3(x, y, z)`                 | From three `Float64` scalars                   |
| `RV3(v::Vector{Float64})`      | From a 3-element Julia vector                  |
| `RV3(v::VEC3)`                 | Reinterpret a `VEC3` as a CRV (no transform)   |
| `RV3()`                        | Identity rotation `(0.0, 0.0, 0.0)`            |
| `RV3(a::RV3, b::RV3)`          | **Compose** two CRVs: R(c) = R(a)·R(b)        |
| `RV3(R::MAT3)`                 | Logarithmic map: rotation matrix → CRV         |

## Cross-type conversion

| Expression    | Result  | Description                               |
|---------------|---------|-------------------------------------------|
| `VEC3(r)`     | `VEC3`  | Reinterpret RV3 as a plain VEC3            |
| `RV3(v)`      | `RV3`   | Reinterpret VEC3 as a CRV                 |

## Scalar operators

Addition and component-wise multiplication between two `RV3` values are
**intentionally not defined** (geometrically meaningless for rotation vectors).

| Expression       | Result  | Description                          |
|------------------|---------|--------------------------------------|
| `r * s`, `s * r` | `RV3`   | Scalar multiplication                |
| `r / s`          | `RV3`   | Scalar division                      |
| `-r`             | `RV3`   | Negation — **reverses** the rotation |

## RV3 geometric operations

| Function           | Returns   | Description                          |
|--------------------|-----------|--------------------------------------|
| `dotp(r1, r2)`     | `Float64` | Dot product  r1 · r2                 |
| `dotp(r, v)`       | `Float64` | Dot product  RV3 · VEC3              |
| `crossp(r1, r2)`   | `RV3`     | Cross product  r1 × r2               |
| `crossp(r, v)`     | `VEC3`    | Cross product  RV3 × VEC3            |
| `norm2(r)`         | `Float64` | Rotation angle ‖r‖ (radians)         |

## RV3 accessors and indexing

| Syntax                    | Description                                    |
|---------------------------|------------------------------------------------|
| `r.x`, `r.y`, `r.z`       | Direct field access                            |
| `r.v`                     | Convert to `Vector{Float64}` `[x, y, z]`       |
| `r[1]`…`r[3]`             | 1-based integer indexing                       |
| `setindex(r, val, i)`     | New `RV3` with component `i` set to `val`      |
| `copy(r)`                 | Deep copy                                      |

```@docs
RV3
norm2(r::RV3)
dotp(r1::RV3, r2::RV3)
crossp(r1::RV3, r2::RV3)
```
