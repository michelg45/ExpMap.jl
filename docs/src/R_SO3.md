```@meta
CurrentModule = ExpMap
```

# R_SO3 — Exponential Map SO(3)

`R_SO3` implements the **exponential map** from the Lie algebra so(3)
(Cartesian Rotation Vectors) to the Lie group SO(3) (3×3 orthogonal rotation
matrices with determinant +1).

## R\_SO3 background

Given a CRV ψ ∈ ℝ³ with rotation angle **θ = ‖ψ‖** and unit axis **k = ψ/θ**,
the rotation matrix is given by the **Rodrigues formula**:

```
R(ψ) = cos θ · I  +  sin θ · k̃  +  (1 − cos θ) · k⊗k
```

where k̃ is the skew-symmetric cross-product matrix of k ([`tilde`](@ref)),
and k⊗k is the outer (dyadic) product ([`outerp`](@ref)).

For small angles (**θ < 1e-8**), the first-order Taylor approximation is used
to avoid numerical division by zero:

```
R(ψ) ≈ I + ψ̃
```

## R\_SO3 methods

### `R_SO3(psi::RV3) → MAT3`

Returns the full 3×3 rotation matrix R(ψ).

### `R_SO3(psi::RV3, a::VEC3) → VEC3`

Rotates the vector `a` by R(ψ) **without forming the matrix**, using the
expanded Rodrigues formula directly on the vector:

```
R(ψ) · a = cos θ · a  +  sin θ · (k × a)  +  (1 − cos θ) · k (k · a)
```

Small-angle approximation:

```
R(ψ) · a ≈ a + ψ × a
```

!!! tip "Efficiency"
    Prefer `R_SO3(psi, a)` over `R_SO3(psi) * a` when only the rotated vector
    is needed — it avoids allocating the full 3×3 matrix.

**Dependencies:** [`tilde`](@ref), [`outerp`](@ref), [`crossp`](@ref), [`dotp`](@ref)

**See also:** [`invR_SO3`](@ref) (logarithmic map), [`NodeFrame`](@ref)

```@docs
R_SO3
```
