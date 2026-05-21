```@meta
CurrentModule = ExpMap
```

# DinvT\_SE3 — Gâteaux Derivative of the Inverse SE(3) Tangent Operator

`DinvT_SE3` computes the 6×6 **matrix of the Gâteaux derivative** of the map
`p ↦ invT_SE3(p)·f` (or its transpose variant) with respect to `p = [u; φ]`,
where `f ∈ ℝ⁶` is a fixed vector.

## DinvT\_SE3 background

The inverse SE(3) tangent operator is:

```
invT_SE3(p) = [ invT    B  ]     B = invT·(−dTu·invT + tilde(Tu))
              [ 0       invT ]
```

Splitting `f = [f_u; f_φ]`, the Gâteaux derivative matrix is:

```
M = [ D_uu  D_uφ ]
    [ D_φu  D_φφ ]
```

## DinvT\_SE3 four cases

| `trp`   | `sign_p` | Computes |
|---------|----------|----------|
| `false` | `"+"`    | `D_p(invT_SE3(p)·f)`           — forward, no transpose |
| `true`  | `"+"`    | `D_p(invT_SE3(p)ᵀ·f)`          — forward, transposed |
| `false` | `"-"`    | `D_p(invT_SE3(p,"-")·f)`       — uses `(invT)ᵀ` base |
| `true`  | `"-"`    | `D_p(invT_SE3(p,"-")ᵀ·f)`      — uses `(invT)ᵀ` base, transposed |

## DinvT\_SE3 finite-difference form

A second method `DinvT_SE3(p::VEC6, f::VEC6)` computes the same matrix by
centred finite differences (step `1e-6`), scaled by `‖u‖` for the translational
components.  Use this for numerical verification:

```julia
p  = VEC6(1.0, -0.5, 0.2,  0.3, -0.1, 0.8)
f  = VEC6(0.1,  0.2, 0.3,  0.4,  0.5, 0.6)
a  = T_SE3_input_data(p)
M1 = DinvT_SE3(a, f)            # analytical
M2 = DinvT_SE3(p, f)            # finite-difference check
```

**Dependencies:** [`T_SE3_data`](@ref), [`tilde`](@ref), [`invT_SE3`](@ref)

**See also:** [`invT_SE3`](@ref), [`DT_SE3`](@ref), [`T_SE3_input_data`](@ref)

```@docs
DinvT_SE3
```
