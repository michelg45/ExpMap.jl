```@meta
CurrentModule = ExpMap
```

# DinvT_SO3 — Gâteaux Derivative of the Inverse Tangent Operator

`DinvT_SO3` computes the 3×3 **matrix of the Gâteaux derivative of the map**
`φ ↦ T⁻¹(φ)·a` (or `φ ↦ (T⁻¹)ᵀ(φ)·a` with option `trp = true`) with respect
to the rotation vector `φ`, where `T⁻¹(φ)` is the inverse SO(3) tangent operator
([`invT_SO3`](@ref)) and `a ∈ ℝ³` is a fixed vector.

## DinvT\_SO3 background

Recall the inverse tangent operator:

```
T⁻¹(φ) = I + ½ φ̃ + γ(x) (φ̂⊗φ̂ − I)
```

where `x = ‖φ‖`, `φ̂ = φ/x`, and `γ(x) = 1 − (x/2) cot(x/2)` (see
[`invT_functions`](@ref)).

The key structural difference from [`DT_SO3`](@ref) is that the linear term
uses `φ` directly (not the normalised `φ̂`).  Its Gâteaux derivative is
therefore the **constant matrix** `−½ã`:

```
d(½ φ×a)/dφ · δφ = ½ (δφ×a) = −½ ã · δφ
```

The `γ`-term is differentiated by the same chain rule as in `DT_SO3`:

```
d(γ·(φ̂⊗φ̂−I)·a)/dφ = γ·DWa·dφ̂ + γ′·Wa⊗φ̂
```

Combining both contributions:

```
d(T⁻¹·a)/dφ = −½·ã + γ·DWa·dφ̂ + γ′·Wa⊗φ̂
```

with:

| Quantity | Definition | Dimension |
|----------|------------|-----------|
| `φ̂`    | `φ/‖φ‖` — unit rotation axis | ℝ³ |
| `dφ̂`   | `(I − φ̂⊗φ̂)/‖φ‖` — Gâteaux derivative of `φ → φ̂` | ℝ³ˣ³ |
| `φ̂`    | gradient of `‖φ‖` (same vector) | ℝ³ |
| `Wa`    | `(φ̂⊗φ̂ − I)·a` — γ-contribution vector | ℝ³ |
| `DWa`   | `φ̂⊗a + (φ̂·a)·I` — Gâteaux derivative of `(φ̂⊗φ̂)·a` w.r.t. `φ̂` | ℝ³ˣ³ |
| `γ′`   | derivative of `γ` w.r.t. `‖φ‖` | ℝ |

## DinvT\_SO3 transpose option

Setting `trp = true` computes instead the Gâteaux derivative of
`φ ↦ (T⁻¹)ᵀ(φ)·a`.  Since the transpose reverses the sign of the
skew-symmetric term:

```
(T⁻¹)ᵀ(φ) = I − ½ φ̃ + γ(x) (φ̂⊗φ̂ − I)
```

only the constant term changes sign:

```
d((T⁻¹)ᵀ·a)/dφ = +½·ã + γ·DWa·dφ̂ + γ′·Wa⊗φ̂
```

Both cases are handled by a single sign parameter `s = trp ? −1 : +1`:

```
d(T̂⁻¹·a)/dφ = −s·½·ã + γ·DWa·dφ̂ + γ′·Wa⊗φ̂
```

## DinvT\_SO3 small-angle limit

For `‖φ‖² < 1e-8`:

```
T⁻¹(φ)  ≈ I + ½φ̃   ⟹   DinvT_SO3(φ, a)           ≈  −½ã   (trp = false)
(T⁻¹)ᵀ  ≈ I − ½φ̃   ⟹   DinvT_SO3(φ, a; trp=true)  ≈  +½ã   (trp = true)
```

This first-order approximation is consistent with the threshold used by
[`invT_SO3`](@ref).

## DinvT\_SO3 role in SO(3) kinematics

`DinvT_SO3` arises in the **tangent stiffness** of SO(3)-parameterised structural
models, in particular wherever strain measures are expressed through `T⁻¹`.
Specifically, the directional derivative of the map `φ ↦ T⁻¹(φ)·u` in a
direction `δφ` reads:

```
δ[T⁻¹(φ)·u] = DinvT_SO3(φ, u) · δφ
```

**Dependencies:** [`invT_functions`](@ref), [`tilde`](@ref), [`outerp`](@ref),
[`dotp`](@ref)

**See also:** [`invT_SO3`](@ref), [`invT_functions`](@ref), [`DT_SO3`](@ref)

```@docs
DinvT_SO3
```
