```@meta
CurrentModule = ExpMap
```

# mat6\_array — MAT6 Collections to Dense Float64 Arrays

`vec_mat6` and `mat_mat6` convert collections of [`MAT6`](@ref) objects into
standard Julia `Array{Float64,2}` matrices by assembling the 6×6 blocks
side-by-side into a single dense block matrix.

## Block layout

### `vec_mat6` — Vector{MAT6} → 6 × 6N

A length-`N` vector is assembled into a `6 × 6N` block-row matrix.
Block `k` occupies columns `(k-1)·6+1 : k·6`:

```
┌──────────┬──────────┬─────┬──────────┐
│  V[1]    │  V[2]    │ … │  V[N]    │   6 rows
│  (6×6)   │  (6×6)   │   │  (6×6)   │
└──────────┴──────────┴─────┴──────────┘
 cols 1:6   cols 7:12       cols 6N-5:6N
```

### `mat_mat6` — Matrix{MAT6} → 6M × 6N

An `M × N` matrix of blocks is assembled into a `6M × 6N` block matrix.
Block `(p, q)` occupies rows `(p-1)·6+1 : p·6` and
columns `(q-1)·6+1 : q·6`:

```
┌──────────┬──────────┬─────┐
│ M[1,1]   │ M[1,2]   │ … │   rows  1:6
│──────────┼──────────┼─────│
│ M[2,1]   │ M[2,2]   │ … │   rows  7:12
│──────────┼──────────┼─────│
│    ⋮     │    ⋮     │ ⋱ │
└──────────┴──────────┴─────┘
 cols 1:6   cols 7:12
```

## Usage notes

The resulting `Array{Float64,2}` can be passed directly to Julia's standard
linear algebra routines (`\`, `eigvals`, …) or written to disk via standard
I/O packages.

**Dependencies:** [`MAT6`](@ref)

**See also:** [`MAT6`](@ref), [`VEC6`](@ref)

```@docs
vec_mat6
mat_mat6
```
