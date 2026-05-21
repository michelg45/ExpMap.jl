```@meta
CurrentModule = ExpMap
```

# MAT3 — Immutable 3×3 Matrix

`MAT3` is an **immutable** 3×3 matrix with `Float64` entries stored in
**row-major order**.  Fields are named `aij` where `i` is the row and `j`
the column (e.g. `a12` → row 1, column 2).

Designed to interoperate with [`VEC3`](@ref): the matrix–vector product
`m * v` returns a `VEC3`.

!!! note "Dependencies"
    `MAT3` requires `using LinearAlgebra` (for `diag`, `inv`, `\`) and
    `using Printf` (for `show`).

## MAT3 construction

| Syntax                                | Description                              |
|---------------------------------------|------------------------------------------|
| `MAT3(a11,…,a33)`                     | From nine `Float64` scalars (row-major)  |
| `MAT3(m::Array{Float64,2})`           | From a standard 3×3 Julia matrix         |
| `MAT3(a::VEC3, b::VEC3, c::VEC3)`     | From three **column** vectors            |
| `MAT3()` / `zero(MAT3)`               | Zero matrix                              |
| `one(MAT3)`                           | Identity matrix I₃                      |

## MAT3 arithmetic operators

| Expression        | Result  | Description                          |
|-------------------|---------|--------------------------------------|
| `m1 + m2`         | `MAT3`  | Element-wise addition                |
| `m1 - m2`         | `MAT3`  | Element-wise subtraction             |
| `-m`              | `MAT3`  | Unary negation                       |
| `m1 * m2`         | `MAT3`  | Matrix–matrix product                |
| `m * v`           | `VEC3`  | Matrix–vector product                |
| `s * m`, `m * s`  | `MAT3`  | Scalar multiplication                |
| `m / s`           | `MAT3`  | Scalar division                      |

## MAT3 linear algebra

| Function / Operator  | Returns  | Description                     |
|----------------------|----------|---------------------------------|
| `transpose(m)`       | `MAT3`   | Transpose mᵀ                    |
| `inv(m)`             | `MAT3`   | Matrix inverse m⁻¹              |
| `m \ v`              | `VEC3`   | Solve linear system m·x = v     |
| `diag(v::VEC3)`      | `MAT3`   | Diagonal matrix from a `VEC3`   |

## MAT3 accessors and indexing

| Syntax                      | Description                                    |
|-----------------------------|------------------------------------------------|
| `m.aij`                     | Direct field access (e.g. `m.a12`)             |
| `m.mat`                     | Convert to `Array{Float64,2}` (3×3)            |
| `m[i, j]`                   | Scalar element at row `i`, column `j`          |
| `m[:, j]`                   | Column `j` as a `VEC3`                         |
| `setindex(m, val, i, j)`    | New `MAT3` with element `(i,j)` set to `val`   |
| `copy(m)`                   | Deep copy                                      |

```@docs
MAT3
```
