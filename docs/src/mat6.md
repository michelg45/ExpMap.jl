```@meta
CurrentModule = ExpMap
```

# MAT6 — Immutable 6×6 Matrix

`MAT6` is an **immutable** 6×6 matrix stored as **six column vectors** of type
[`VEC6`](@ref).  This column-major layout makes matrix–vector products
allocation-free and cache-friendly.

In rigid-body mechanics a `MAT6` typically represents a **spatial operator**
such as a stiffness, compliance, or inertia matrix acting on wrenches or twists.

!!! note "Dependencies"
    `MAT6` requires `using LinearAlgebra` and `using Printf`.
    `inv` and `\` delegate to Julia's built-in dense solver; the matrix must
    be non-singular.

## MAT6 construction

| Syntax                                  | Description                                      |
|-----------------------------------------|--------------------------------------------------|
| `MAT6(c1,…,c6::VEC6)`                  | From six column `VEC6` vectors                   |
| `MAT6(v::Vector{Float64})`              | From a 36-element flat vector (column-major)     |
| `MAT6(A::Array{Float64,2})`             | From a standard 6×6 Julia matrix                 |
| `MAT6(A, B, C, D ::MAT3)`              | Block form `[A B; C D]` from four 3×3 blocks     |
| `MAT6()` / `zero(MAT6)`                 | Zero matrix                                      |
| `one(MAT6)`                             | Identity matrix I₆                              |

The block constructor layout is:

```
┌         ┐
│  A  │ B │   rows 1–3,  cols 1–3 / 4–6
│─────┼───│
│  C  │ D │   rows 4–6
└         ┘
```

## MAT6 arithmetic operators

| Expression        | Result  | Description               |
|-------------------|---------|---------------------------|
| `m1 + m2`         | `MAT6`  | Column-wise addition      |
| `m1 - m2`         | `MAT6`  | Column-wise subtraction   |
| `m * v`           | `VEC6`  | Matrix–vector product     |
| `m1 * m2`         | `MAT6`  | Matrix–matrix product     |
| `s * m`, `m * s`  | `MAT6`  | Scalar multiplication     |

## MAT6 linear algebra

| Function / Operator | Returns  | Description                           |
|---------------------|----------|---------------------------------------|
| `transpose(m)`      | `MAT6`   | Transpose mᵀ                          |
| `inv(m)`            | `MAT6`   | Matrix inverse m⁻¹                    |
| `m \ v`             | `VEC6`   | Solve m·x = v                         |
| `m \ n`             | `MAT6`   | Solve m·X = n (column by column)      |
| `diag(v::VEC6)`     | `MAT6`   | Diagonal matrix from a `VEC6`         |

## MAT6 accessors and indexing

| Syntax                    | Description                                        |
|---------------------------|----------------------------------------------------|
| `m.v1` … `m.v6`           | Column `j` as a `VEC6`                             |
| `m.mat`                   | Convert to `Array{Float64,2}` (6×6)               |
| `m[j]`                    | Column `j` as `VEC6` (1-based)                     |
| `m[i, j]`                 | Scalar element at row `i`, column `j`              |
| `setindex(m, val, i, j)`  | New `MAT6` with element `(i,j)` set to `val`       |

```@docs
MAT6
```
