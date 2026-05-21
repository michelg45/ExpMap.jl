```@meta
CurrentModule = ExpMap
```

# invT_SO3 — Inverse Tangent Operator T⁻¹(φ)

`invT_SO3` computes the **inverse tangent operator** T⁻¹(φ) of SO(3),
or its action on a vector `a`, without forming the full matrix.

## invT\_SO3 background

The inverse tangent operator is the left-trivialized inverse of the tangent
map of the exponential on SO(3).  For a rotation vector φ ∈ ℝ³ with norm
x = ‖φ‖, it reads:

```
T⁻¹(φ) = I + ½ φ̃ + γ(x) (φ̂ ⊗ φ̂ − I)
```

where:
- **φ̂ = φ / x** — unit rotation axis
- **φ̃** — skew-symmetric matrix of φ ([`tilde`](@ref))
- **γ(x) = 1 − (x/2) cot(x/2)** — scalar coefficient ([`invT_functions`](@ref))

The relation `T⁻¹(φ) T(φ) = I` holds for all φ.

## invT\_SO3 small-angle limit

For ‖φ‖ < 1e-4 the first-order approximation is returned:

```
T⁻¹(φ) ≈ I + ½ φ̃
```

The vector form correspondingly gives `T⁻¹(φ) a ≈ a + ½ (φ × a)`.

## invT\_SO3 methods

### `invT_SO3(φ::RV3) → MAT3`

Returns the full 3×3 inverse tangent matrix T⁻¹(φ).

### `invT_SO3(φ::RV3, a::VEC3) → VEC3`

Returns T⁻¹(φ) a **without forming the full matrix**:

```
T⁻¹(φ) a = a + ½ (φ × a) + (γ/x²) (φ (φ · a) − x² a)
```

**Dependencies:** [`invT_functions`](@ref), [`tilde`](@ref), [`outerp`](@ref), [`crossp`](@ref), [`dotp`](@ref)

**See also:** [`invT_functions`](@ref), [`T_functions`](@ref)

```@docs
invT_SO3
```
