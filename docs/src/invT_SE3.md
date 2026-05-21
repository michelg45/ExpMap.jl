```@meta
CurrentModule = ExpMap
```

# invT\_SE3 — Inverse SE(3) Tangent Operator

`invT_SE3` constructs the 6×6 inverse SE(3) tangent operator `T_SE3(p)⁻¹`
from either a `VEC6` parameter vector or a precomputed [`T_SE3_data`](@ref).

## invT\_SE3 background

For `p = [u; φ]`, exploiting the block upper-triangular structure of
[`T_SE3`](@ref):

```
invT_SE3(p) = [ T⁻¹(φ)    T⁻¹(φ)·(−dTu·T⁻¹(φ) + tilde(Tu)) ]
              [ 0₃ₓ₃       T⁻¹(φ)                             ]
```

with `invT = T⁻¹(φ)`, `Tu = T(φ)·u`, `dTu = DT_SO3(φ, u)`.

## invT\_SE3 transpose variant

With `trp = true` (requires [`T_SE3_data`](@ref) input), the substitutions are:

```
invT         →   (T⁻¹)ᵀ(φ)
Tu           →  −(T⁻¹)ᵀ(φ)·u
dTu[:,i]     →   (∂T/∂φᵢ)ᵀ·u
```

## invT\_SE3 identity check

```julia
p  = VEC6(1.0, -0.5, 0.2,  0.3, -0.1, 0.8)
a  = T_SE3_input_data(p)
T_SE3(a) * invT_SE3(a)   # should equal I₆  (identity MAT6)
```

**Dependencies:** [`invT_SE3_input`](@ref)

**See also:** [`T_SE3`](@ref), [`DinvT_SE3`](@ref)

```@docs
invT_SE3
```
