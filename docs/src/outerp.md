```@meta
CurrentModule = ExpMap
```

# outerp — Outer (Dyadic) Product

`outerp(u, v)` computes the **outer product** (also called the dyadic or
tensor product) of two 3D vectors, returning a 3×3 matrix.

## outerp definition

Given **u, v ∈ ℝ³**, the outer product is the rank-1 matrix:

```
(u ⊗ v)ᵢⱼ = uᵢ · vⱼ

┌ u.x·v.x   u.x·v.y   u.x·v.z ┐
│ u.y·v.x   u.y·v.y   u.y·v.z │
└ u.z·v.x   u.z·v.y   u.z·v.z ┘
```

!!! note
    `outerp(u, v) ≠ outerp(v, u)` in general.
    The result is symmetric only when `u = v`, giving a rank-1 projection matrix.

## outerp role in SO(3)

`outerp` appears in the **Rodrigues rotation formula** ([`R_SO3`](@ref)):

```
R(ψ) = cos θ · I  +  sin θ · k̃  +  (1 − cos θ) · outerp(k, k)
```

where k = ψ/‖ψ‖ is the unit rotation axis.

## outerp methods

| Signature                         | Description                                       |
|-----------------------------------|---------------------------------------------------|
| `outerp(v1::VEC3, v2::VEC3) → MAT3` | Outer product of two `VEC3`                    |
| `outerp(v1::RV3,  v2::RV3)  → MAT3` | Same for `RV3` (avoids explicit conversion)    |

**See also:** [`tilde`](@ref), [`R_SO3`](@ref)

```@docs
outerp
```
