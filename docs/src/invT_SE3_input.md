```@meta
CurrentModule = ExpMap
```

# invT\_SE3\_input — Auxiliary quantities for the inverse SE(3) tangent operator

`invT_SE3_input` computes the three SO(3) building blocks needed to assemble
`T_SE3(p)⁻¹`, sharing all intermediate calculations.

## invT\_SE3\_input background

Given `p = [u; φ]`, the inverse SE(3) tangent operator reads:

```
invT_SE3(p) = [ T⁻¹(φ)    T⁻¹(φ)·(−dTu·T⁻¹(φ) + tilde(Tu)) ]
              [ 0₃ₓ₃       T⁻¹(φ)                             ]
```

The three outputs of `invT_SE3_input` are:

| Output  | Expression      | Rôle |
|---------|-----------------|------|
| `invT`  | `T⁻¹(φ)`        | inverse SO(3) tangent operator |
| `Tu`    | `T(φ)·u`        | tangent image of translation |
| `dTu`   | `DT_SO3(φ, u)`  | Gâteaux derivative of `φ ↦ T(φ)·u` |

Note that the building blocks for the *inverse* SE(3) operator involve `T(φ)·u`
and `DT_SO3(φ,u)` — not their `T⁻¹` counterparts — because the block formula
couples `T⁻¹` with the forward-tangent quantities.

## invT\_SE3\_input small-angle limit

For `‖φ‖² < 1e-8`:

```
invT ≈ I + ½φ̃,   Tu ≈ u − ½(φ×u),   dTu ≈ ½ũ
```

**Dependencies:** [`T_functions`](@ref), [`invT_functions`](@ref),
[`tilde`](@ref), [`outerp`](@ref)

**See also:** [`invT_SE3`](@ref), [`T_SE3_input`](@ref)

```@docs
invT_SE3_input
```
