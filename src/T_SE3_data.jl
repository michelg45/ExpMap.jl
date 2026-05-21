"""
    T_SE3_data

Precomputed data for the SE(3) tangent operator and its directional derivative.
Populated by [`T_SE3_input_data`](@ref); consumed by [`T_SE3`](@ref),
[`invT_SE3`](@ref), [`DT_SE3`](@ref), and [`DinvT_SE3`](@ref).

| Field     | Type              | Content |
|-----------|-------------------|---------|
| `p`       | `VEC6`            | input parameter vector `[u; φ]` |
| `invT`    | `MAT3`            | `T⁻¹(φ)` — inverse SO(3) tangent operator |
| `T`       | `MAT3`            | `T(φ)`   — SO(3) tangent operator |
| `Tu`      | `VEC3`            | `T(φ)·u` |
| `dTu`     | `MAT3`            | Gâteaux derivative matrix of `φ ↦ T(φ)·u` |
| `dinvT`   | `Vector{MAT3}[3]` | `dinvT[i] = ∂T⁻¹/∂φᵢ` |
| `dT`      | `Vector{MAT3}[3]` | `dT[i]   = ∂T/∂φᵢ` |
| `d²Tu`    | `Vector{MAT3}[3]` | `d²Tu[i][:,j]  = (∂²T/∂φᵢ∂φⱼ)·u` |
| `d²Ttu`   | `Vector{MAT3}[3]` | `d²Ttu[i][:,j] = (∂²Tᵀ/∂φᵢ∂φⱼ)·u` |
"""
struct T_SE3_data
    p::VEC6
    invT::MAT3
    T::MAT3
    Tu::VEC3
    dTu::MAT3
    dinvT::Vector{MAT3}
    dT::Vector{MAT3}
    d²Tu::Vector{MAT3}
    d²Ttu::Vector{MAT3}
end