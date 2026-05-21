```@meta
CurrentModule = ExpMap
```

# T_SO3 — Tangent Operator T(φ)

`T_SO3` computes the **tangent operator** T(φ) of SO(3), or its action on a
vector `a`, without forming the full matrix.

## T\_SO3 background

The tangent operator (also called the left-trivialized tangent map) arises when
linearising the exponential map on SO(3).  For a rotation vector φ ∈ ℝ³ with
norm x = ‖φ‖, it reads:

```
T(φ) = I + α(x) φ̂̃ + β(x) (φ̂ ⊗ φ̂ − I)
```

where:
- **φ̂ = φ / x** — unit rotation axis
- **φ̂̃** — skew-symmetric matrix of φ̂ ([`tilde`](@ref))
- **α(x) = (cos x − 1) / x** — scalar coefficient ([`T_functions`](@ref))
- **β(x) = 1 − sin(x) / x** — scalar coefficient ([`T_functions`](@ref))

## T\_SO3 small-angle limit

For ‖φ‖ < 1e-4 the first-order approximation is returned to avoid
numerical cancellation:

```
T(φ) ≈ I − ½ φ̃
```

The vector form correspondingly gives `T(φ) a ≈ a − ½ (φ × a)`.

## T\_SO3 methods

### `T_SO3(φ::RV3) → MAT3`

Returns the full 3×3 tangent matrix T(φ).

### `T_SO3(φ::RV3, a::VEC3) → VEC3`

Returns T(φ) a **without forming the full matrix**:

```
T(φ) a = a + α(x) (φ̂ × a) + β(x) (φ̂ (φ̂ · a) − a)
```

!!! tip "Efficiency"
    Prefer `T_SO3(φ, a)` over `T_SO3(φ) * a` when only the image vector
    is needed — it avoids allocating the full 3×3 matrix.

**Dependencies:** [`T_functions`](@ref), [`tilde`](@ref), [`outerp`](@ref), [`crossp`](@ref), [`dotp`](@ref)

**See also:** [`invT_SO3`](@ref), [`T_functions`](@ref), [`R_SO3`](@ref)

```@docs
T_SO3
```
