```@meta
CurrentModule = ExpMap
```

# shape\_functions\_1D — 1D Lagrange shape functions

`shape_functions_1D` evaluates the 1D Lagrange interpolation polynomials and
their derivatives at a single point `ξ` on the reference segment [−1, +1].

## shape\_functions\_1D elements

| `Nnodes` | Element | Node locations |
|----------|---------|----------------|
| 2 | linear L2     | ξ = −1, +1 |
| 3 | quadratic L3  | ξ = −1, 0, +1 |
| 4 | cubic L4      | ξ = −1, −1/3, +1/3, +1 |

The L4 shape functions involve the factor 9/16 = 0.5625 arising from the
Lagrange product formula evaluated at nodes ±1 and ±1/3.

## shape\_functions\_1D formulas

**L2** (`Nnodes = 2`):
```
N₁ = ½(1−ξ),   N₂ = ½(1+ξ)
```

**L3** (`Nnodes = 3`):
```
N₁ = ½ξ(ξ−1),   N₂ = 1−ξ²,   N₃ = ½ξ(1+ξ)
```

## shape\_functions\_1D partition of unity

For all elements: `∑ F[i] = 1` and `∑ DF[i] = 0` at every `ξ`.

## shape\_functions\_1D example

```julia
# Interpolate a function known at 3 nodes
x, w = gauss_points(1, 3)         # 3 Gauss points
for i = 1:3
    F, DF = shape_functions_1D(3, x[i])   # shape functions at Gauss point i
    # F[j] is the weight of node j at integration point i
end
```

**See also:** [`gauss_points`](@ref), [`shape_functions_2D`](@ref)

```@docs
shape_functions_1D
```
