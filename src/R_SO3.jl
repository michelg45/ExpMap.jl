# ==============================================================================
# R_SO3.jl — SO(3) Rotation Matrix and Rotation-of-Vector Functions
# ==============================================================================
#
# Provides two methods of `R_SO3` that implement the exponential map from the
# Lie algebra so(3) (Cartesian Rotation Vectors) to the Lie group SO(3)
# (3×3 orthogonal rotation matrices with determinant +1).
#
# Given a Cartesian Rotation Vector (CRV)  ψ ∈ ℝ³:
#   • θ = ‖ψ‖   is the rotation angle  (radians)
#   • k = ψ/θ   is the unit rotation axis  (when θ > 0)
#
# The rotation matrix is given by the Rodrigues formula:
#
#   R(ψ) = cos θ · I  +  sin θ · k̃  +  (1 − cos θ) · k⊗k
#
# where  k̃  is the skew-symmetric cross-product matrix of  k,  and  k⊗k
# is the outer (dyadic) product.
#
# For small angles (θ < 1e-8), a first-order approximation is used to avoid
# numerical division by zero:
#
#   R(ψ) ≈ I + ψ̃        (first-order Taylor expansion)
#
# METHODS
#   R_SO3(psi::RV3)           → MAT3   full 3×3 rotation matrix R(ψ)
#   R_SO3(psi::RV3, a::VEC3) → VEC3   rotate vector a by R(ψ)  (avoids
#                                       building the full matrix)
#
# RODRIGUES FORMULA — EXPANDED FORM  (second method)
#   R(ψ) · a = cos θ · a  +  sin θ · (k × a)  +  (1 − cos θ) · k (k · a)
#
#   Small-angle approximation:
#   R(ψ) · a ≈ a + ψ × a
#
# EXTERNAL DEPENDENCIES
#   tilde(v::VEC3)            → MAT3   skew-symmetric matrix of v  (ṽ · w = v × w)
#   outerp(u::VEC3, v::VEC3) → MAT3   outer (dyadic) product  u ⊗ v
#   crossp(u::VEC3, v::VEC3) → VEC3   cross product  u × v
#   dotp(u::RV3,  v::RV3)   → Float64 dot product
#   dotp(u::VEC3, v::VEC3)  → Float64 dot product
#   one(MAT3)                → MAT3   identity matrix  I₃
#   VEC3, RV3, MAT3          see VEC3.jl, RV3.jl, MAT3.jl
#
# NOTES
#   • The small-angle threshold (1e-8 rad) safely avoids division by zero
#     while keeping the first-order truncation error below machine precision
#     for Float64.
#   • The second method `R_SO3(psi, a)` is more efficient than
#     `R_SO3(psi) * a` when only the rotated vector is needed, as it
#     requires no matrix allocation.
# ==============================================================================


"""
    R_SO3(psi::RV3) → MAT3

Return the 3×3 rotation matrix `R(ψ)` corresponding to the Cartesian
Rotation Vector `psi`, via the Rodrigues formula:

    R(ψ) = cos θ · I  +  sin θ · k̃  +  (1 − cos θ) · k⊗k

where `θ = ‖ψ‖` is the rotation angle and `k = ψ/θ` is the unit axis.

For small angles (`θ < 1e-8`), the first-order approximation `R ≈ I + ψ̃`
is used to avoid numerical instability.
"""
function R_SO3(psi::RV3)
    theta = sqrt(dotp(psi, psi))   # rotation angle θ = ‖ψ‖

    if theta < 1e-8
        # First-order Taylor expansion: R(ψ) ≈ I + ψ̃
        return one(MAT3) + tilde(psi)
    else
        k = psi / theta   # unit rotation axis
        # Full Rodrigues formula: R = cos θ·I + sin θ·k̃ + (1−cos θ)·k⊗k
        return cos(theta) * one(MAT3) + sin(theta) * tilde(k) + (1 - cos(theta)) * outerp(k, k)
    end
end


"""
    R_SO3(psi::RV3, a::VEC3) → VEC3

Rotate the vector `a` by the SO(3) rotation defined by the CRV `psi`,
using the expanded Rodrigues formula directly on the vector:

    R(ψ) · a = cos θ · a  +  sin θ · (k × a)  +  (1 − cos θ) · k (k · a)

where `θ = ‖ψ‖` and `k = ψ/θ`.

For small angles (`θ < 1e-8`), the first-order approximation is used:

    R(ψ) · a ≈ a + ψ × a

This method is more efficient than `R_SO3(psi) * a` as it avoids
constructing the full 3×3 rotation matrix.
"""
function R_SO3(psi::RV3, a::VEC3)
    theta = sqrt(dotp(psi, psi))   # rotation angle θ = ‖ψ‖

    if theta < 1e-8
        # First-order approximation: R(ψ)·a ≈ a + ψ × a
        return a + crossp(VEC3(psi), a)
    else
        k = VEC3(psi[1] / theta, psi[2] / theta, psi[3] / theta)   # unit axis
        # Expanded Rodrigues: cos θ·a + sin θ·(k×a) + (1−cos θ)·k(k·a)
        return cos(theta) * a + sin(theta) * crossp(k, a) + (1 - cos(theta)) * k * dotp(k, a)
    end
end

"""
    R_SO3(phi::Float64, axe::Int) → MAT3

Return the 3×3 rotation matrix for a rotation of angle `phi` (radians)
about coordinate axis `axe` (1 = x-axis, 2 = y-axis, 3 = z-axis).

Equivalent to `R_SO3(RV3(phi·eₐₓₑ))` where `eₐₓₑ` is the unit basis vector
along axis `axe`.

## Example
```julia
julia> R_SO3(π/2, 3)    # 90° rotation about the z-axis
```

**See also:** [`R_SO3(psi::RV3)`](@ref), [`euler_to_RV3`](@ref)
"""
function R_SO3(phi::Float64, axe::Int)
    n      = zeros(3)
    n[axe] = 1.0
    return R_SO3(RV3(phi * n))
end