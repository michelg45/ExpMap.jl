```@meta
CurrentModule = ExpMap
```

# shape\_functions\_2D — 2D Lagrange shape functions

`shape_functions_2D` evaluates 2D Lagrange interpolation polynomials and their
partial derivatives at a point `ξ = [ξ₁, ξ₂]` on a 2D reference domain.

## shape\_functions\_2D elements

Two reference domains are used depending on the element type.

### Reference triangle — vertices (0,0), (1,0), (0,1)

| `Nnodes` | Element | Node layout |
|----------|---------|-------------|
| 3 | linear triangle T3   | 3 vertices |
| 6 | quadratic triangle T6 | 3 vertices + 3 mid-side nodes at (½,0),(½,½),(0,½) |

### Reference square — [−1,+1]²

| `Nnodes` | Element | Node layout |
|----------|---------|-------------|
| 4 | bilinear quad Q1          | 4 corners counter-clockwise |
| 5 | enriched Q1 + bubble node | 4 corners + central node at (0,0) |
| 8 | serendipity Q2s           | 4 corners + 4 mid-side nodes |
| 9 | biquadratic Q2 (Lagrange) | 4 corners + 4 mid-side + 1 centre |

Corners are numbered counter-clockwise starting at (−1,−1).

## shape\_functions\_2D output

`DF` is an `Nnodes × 2` matrix: `DF[i,1] = ∂Fᵢ/∂ξ₁`, `DF[i,2] = ∂Fᵢ/∂ξ₂`.

For an isoparametric mapping `x(ξ)`, the Jacobian is `J = DF' * x_nodes`
(2×2 matrix).

## shape\_functions\_2D partition of unity

For all elements: `∑ F[i] = 1` and `∑ DF[i,:] = [0 0]` at every `ξ`.

## shape\_functions\_2D example

```julia
# Integrate over a bilinear quad element
x, w = gauss_points(2, 4)         # 2×2 Gauss rule
for i = 1:4
    F, DF = shape_functions_2D(4, x[i])
    J  = DF' * x_nodes             # 2×2 Jacobian (x_nodes: 4×2 node coords)
    dV = w[i] * abs(det(J))        # integration weight
end
```

**See also:** [`gauss_points`](@ref), [`shape_functions_1D`](@ref)

```@docs
shape_functions_2D
```
