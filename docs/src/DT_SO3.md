```@meta
CurrentModule = ExpMap
```

# DT_SO3 — Gâteaux Derivative of the Tangent Operator

`DT_SO3` computes the 3×3 **matrix of the Gâteaux derivative of the map**
`φ ↦ T(φ)·a` (or `φ ↦ Tᵀ(φ)·a` with option `trp = true`) with respect
to the rotation vector `φ`, where `T(φ)` is the SO(3) tangent operator
([`T_SO3`](@ref)) and `a ∈ ℝ³` is a fixed vector.

## DT\_SO3 background

Recall the tangent operator:

```
T(φ) = I + α(x) φ̂̃ + β(x) (φ̂⊗φ̂ − I)
```

where `x = ‖φ‖`, `φ̂ = φ/x`, and `α, β` are scalar coefficients (see
[`T_functions`](@ref)).  Differentiating `T(φ)·a` with respect to `φ` via
the chain rule gives:

```
d(T·a)/dφ = [−α·ã + β·DWa] · dφ̂  +  α′·(φ̂×a)⊗φ̂  +  β′·Wa⊗φ̂
```

with:

| Quantity | Definition | Dimension |
|----------|------------|-----------|
| `φ̂`    | `φ/‖φ‖` — unit rotation axis | ℝ³ |
| `dφ̂`   | `(I − φ̂⊗φ̂)/‖φ‖` — Gâteaux derivative of `φ → φ̂` | ℝ³ˣ³ |
| `φ̂`    | gradient of `‖φ‖` (same vector) | ℝ³ |
| `Wa`    | `(φ̂⊗φ̂ − I)·a` — β-contribution vector | ℝ³ |
| `DWa`   | `φ̂⊗a + (φ̂·a)·I` — Gâteaux derivative of `(φ̂⊗φ̂)·a` w.r.t. `φ̂` | ℝ³ˣ³ |
| `α′, β′`| derivatives of `α, β` w.r.t. `‖φ‖` | ℝ |

## DT\_SO3 transpose option

Setting `trp = true` computes instead the Gâteaux derivative of
`φ ↦ Tᵀ(φ)·a`.  Since the transpose reverses the sign of the skew-symmetric
term:

```
Tᵀ(φ) = I − α(x) φ̃ + β(x) (φ̂⊗φ̂ − I)
```

the only change in the formula is a sign flip on the `α`-dependent terms:

```
d(Tᵀ·a)/dφ = [+α·ã + β·DWa] · dφ̂  −  α′·(φ̂×a)⊗φ̂  +  β′·Wa⊗φ̂
```

Both cases are handled by a single sign parameter `s = trp ? −1 : +1`:

```
(−s·α·ã + β·DWa) · dφ̂  +  s·α′·(φ̂×a)⊗φ̂  +  β′·Wa⊗φ̂
```

## DT\_SO3 small-angle limit

For `‖φ‖² < 1e-8`:

```
T(φ)  ≈ I − ½φ̃   ⟹   DT_SO3(φ, a)          ≈  +½ã   (trp = false)
Tᵀ(φ) ≈ I + ½φ̃   ⟹   DT_SO3(φ, a; trp=true) ≈  −½ã   (trp = true)
```

This first-order approximation is consistent with the threshold used by
[`T_SO3`](@ref).

## DT\_SO3 role in SO(3) kinematics

`DT_SO3` arises in the **tangent stiffness** of SO(3)-parameterised structural
models.  Specifically, when computing the Gâteaux derivative of the map
`φ ↦ T(φ)·u` (with `u` a fixed vector) in a direction `δφ`, the result is:

```
δ[T(φ)·u] = DT_SO3(φ, u) · δφ
```

This matrix is needed in second-order linearisations of the rotation kinematics,
for instance in the linearisation of the tangent operator itself with respect
to the rotation parameter.

**Dependencies:** [`T_functions`](@ref), [`tilde`](@ref), [`outerp`](@ref),
[`crossp`](@ref), [`dotp`](@ref)

**See also:** [`T_SO3`](@ref), [`T_functions`](@ref), [`bracket`](@ref)

```@docs
DT_SO3
```
