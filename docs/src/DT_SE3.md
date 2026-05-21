```@meta
CurrentModule = ExpMap
```

# DT\_SE3 — Gâteaux Derivative of the SE(3) Tangent Operator

`DT_SE3` computes the 6×6 **matrix of the Gâteaux derivative** of the map
`p ↦ T_SE3(p)·f` (or its transpose variant) with respect to the full parameter
vector `p = [u; φ]`, where `f ∈ ℝ⁶` is a fixed vector.

## DT\_SE3 background

The SE(3) tangent operator is block upper-triangular:

```
T_SE3(p) = [ T      B  ]     B = dTu − tilde(Tu)·T
           [ 0      T  ]
```

Splitting `f = [f_u; f_φ]`, the action `T_SE3(p)·f` has:

- **lower block**: `T·f_φ`  — depends only on `φ`
- **upper block**: `T·f_u + B·f_φ`  — depends on both `u` and `φ`

The Gâteaux derivative is the 6×6 matrix `M` such that
`δ(T_SE3(p)·f) = M·δp`.  Written in `[u; φ]` blocks:

```
M = [ D_uu  D_uφ ]
    [ D_φu  D_φφ ]
```

## DT\_SE3 four cases

| `trp`   | `sign_p` | Computes |
|---------|----------|----------|
| `false` | `"+"`    | `D_p(T_SE3(p)·f)`    — forward, no transpose |
| `true`  | `"+"`    | `D_p(T_SE3(p)ᵀ·f)`   — forward, transposed |
| `false` | `"-"`    | `D_p(T_SE3(p,"-")·f)` — uses `Tᵀ` base |
| `true`  | `"-"`    | `D_p(T_SE3(p,"-")ᵀ·f)` — uses `Tᵀ` base, transposed |

## DT\_SE3 block structure (case `!trp, "+"`)

```
D_φu = 0,                      D_φφ = pr3d(dT, f_φ)   [= DT_SO3(φ, f_φ)]
D_uu = D_φφ                    (tilde term omitted, see source header)
D_uφ = pr3d(dT, f_u) + pr3d(d²Tu, f_φ) − tilde(Tu)·D_φφ + tilde(T·f_φ)·dTu
```

where `pr3d(M, v)` denotes the MAT3 whose `j`-th column is `M[j]·v`.

## DT\_SE3 example

```julia
p = VEC6(1.0, 0.0, 0.0,  0.3, -0.1, 0.5)
f = VEC6(0.0, 0.0, 1.0,  1.0,  0.0, 0.0)
a = T_SE3_input_data(p)
M = DT_SE3(a, f)                    # 6×6 Gâteaux derivative matrix
```

**Dependencies:** [`T_SE3_data`](@ref), [`tilde`](@ref)

**See also:** [`T_SE3`](@ref), [`DinvT_SE3`](@ref), [`T_SE3_input_data`](@ref)

```@docs
DT_SE3
```
