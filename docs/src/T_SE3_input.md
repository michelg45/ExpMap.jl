```@meta
CurrentModule = ExpMap
```

# T\_SE3\_input — Auxiliary quantities for the SE(3) tangent operator

`T_SE3_input` computes the three SO(3) building blocks needed to assemble the
SE(3) tangent operator [`T_SE3`](@ref), sharing all intermediate calculations
so that `φ̂`, `α`, `β`, `α′`, `β′`, and `φ̂⊗φ̂` are evaluated only once.

## T\_SE3\_input background

Given the parameter vector `p = [u; φ]` with translation `u ∈ ℝ³` and rotation
vector `φ ∈ ℝ³`, the three outputs are:

```
T(φ)   = I + α(x) tilde(φ̂) + β(x) W           SO(3) tangent operator
Tu     = u + α(x) (φ̂×u) + β(x) W·u             T(φ)·u
dTu    = (−α ũ + β DWu) · dφ̂ + α′ (φ̂×u)⊗φ̂ + β′ Wu⊗φ̂
```

with `x = ‖φ‖`, `φ̂ = φ/x`, `W = φ̂⊗φ̂ − I`, `dφ̂ = (I − φ̂⊗φ̂)/x`.
See [`T_SO3`](@ref) and [`DT_SO3`](@ref) for the scalar coefficients.

## T\_SE3\_input small-angle limit

For `‖φ‖² < 1e-8`:

```
T   ≈ I − ½φ̃
Tu  ≈ u − ½(φ×u)
dTu ≈ ½ũ
```

## T\_SE3\_input role

This function is the lightweight entry point for assembling `T_SE3(p)` from a
`VEC6` without pre-allocating a full [`T_SE3_data`](@ref) struct.  Use
[`T_SE3_input_data`](@ref) instead when the directional derivative
[`DT_SE3`](@ref) or [`DinvT_SE3`](@ref) is also needed.

**Dependencies:** [`T_functions`](@ref), [`tilde`](@ref), [`outerp`](@ref)

**See also:** [`T_SE3`](@ref), [`T_SE3_input_data`](@ref)

```@docs
T_SE3_input
```
