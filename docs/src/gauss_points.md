```@meta
CurrentModule = ExpMap
```

# gauss\_points — Gauss–Legendre quadrature rules

`gauss_points` returns the abscissae and weights for Gauss–Legendre quadrature
on standard reference domains in 1D and 2D.

## gauss\_points reference domains

### 1D — segment [−1, +1]

| `npoints` | Element | Exact for polynomials up to |
|-----------|---------|-----------------------------|
| 1         | midpoint | degree 1 |
| 2         | 2-point  | degree 3 |
| 3         | 3-point  | degree 5 |
| 4         | 4-point  | degree 7 |

Weights sum to 2 (length of [−1,+1]).

### 2D — square [−1,+1]² (tensor-product rules)

| `npoints` | Rule | Notes |
|-----------|------|-------|
| 1  | 1-point  | centroid |
| 4  | 2×2 Gauss–Legendre | |
| 9  | 3×3 Gauss–Legendre | |

Weights sum to 4 (area of [−1,+1]²).

### 2D — reference triangle, vertices (0,0),(1,0),(0,1)

| `npoints` | Rule | Notes |
|-----------|------|-------|
| 7  | Dunavant degree-5 | weights scaled by ½, sum to ½ |

## gauss\_points output format

For `ndim = 1`, `x` is a `Vector{Float64}` and `w` a `Vector{Float64}`.

For `ndim = 2`, `x` is a `Vector{Vector{Float64}}` (each inner vector is a
2-component coordinate) and `w` a `Vector{Float64}`.

## gauss\_points example

```julia
# 1D — 3-point rule
x, w = gauss_points(1, 3)   # x[i] ∈ [−1,1], sum(w) = 2

# 2D — 2×2 rule on the square
x, w = gauss_points(2, 4)   # x[i] = [ξ₁,ξ₂], sum(w) = 4

# Integrate f(ξ) = ξ² on [−1,+1]:
x, w = gauss_points(1, 2)
I = sum(w[i] * x[i]^2 for i = 1:2)   # ≈ 2/3
```

**See also:** [`shape_functions_1D`](@ref), [`shape_functions_2D`](@ref)

```@docs
gauss_points
```
