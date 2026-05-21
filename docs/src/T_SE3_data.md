```@meta
CurrentModule = ExpMap
```

# T\_SE3\_data — Data structure for the SE(3) tangent operator

`T_SE3_data` is an **immutable struct** that bundles all SO(3) quantities
needed to assemble the SE(3) tangent operator and its Gâteaux derivative.
It is produced by [`T_SE3_input_data`](@ref) and consumed by [`T_SE3`](@ref),
[`invT_SE3`](@ref), [`DT_SE3`](@ref), and [`DinvT_SE3`](@ref).

## T\_SE3\_data fields

| Field    | Type              | Meaning |
|----------|-------------------|---------|
| `p`      | `VEC6`            | Input parameter vector `[u; φ]` |
| `invT`   | `MAT3`            | `T⁻¹(φ)` — inverse SO(3) tangent operator |
| `T`      | `MAT3`            | `T(φ)`   — SO(3) tangent operator |
| `Tu`     | `VEC3`            | `T(φ)·u` |
| `dTu`    | `MAT3`            | Gâteaux derivative matrix of `φ ↦ T(φ)·u` |
| `dinvT`  | `Vector{MAT3}[3]` | `dinvT[i] = ∂T⁻¹/∂φᵢ` |
| `dT`     | `Vector{MAT3}[3]` | `dT[i]   = ∂T/∂φᵢ` |
| `d²Tu`   | `Vector{MAT3}[3]` | `d²Tu[i][:,j]  = (∂²T/∂φᵢ∂φⱼ)·u` |
| `d²Ttu`  | `Vector{MAT3}[3]` | `d²Ttu[i][:,j] = (∂²Tᵀ/∂φᵢ∂φⱼ)·u` |

## T\_SE3\_data usage pattern

```julia
# Populate once
a = T_SE3_input_data(p)

# Assemble operators from precomputed data
M  = T_SE3(a)                             # forward operator
Mi = invT_SE3(a)                          # inverse operator
DM = DT_SE3(a, f)                         # directional derivative of T_SE3·f
DMi = DinvT_SE3(a, f)                     # directional derivative of invT_SE3·f
```

**See also:** [`T_SE3_input_data`](@ref), [`T_SE3`](@ref), [`invT_SE3`](@ref)

```@docs
T_SE3_data
```
