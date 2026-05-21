# ==============================================================================
# tilde.jl — Skew-Symmetric (Cross-Product) Matrix of a 3D Vector
# ==============================================================================
#
# Provides two methods of `tilde` that construct the 3×3 skew-symmetric matrix
# associated with a 3D vector, also known as the hat map or cross-product matrix.
#
# DEFINITION
#   Given a vector  v = (vx, vy, vz)ᵀ,  the skew-symmetric matrix is:
#
#       ṽ  =  tilde(v)  =  ┌  0   -vz   vy ┐
#                           │  vz   0   -vx │
#                           └ -vy   vx   0  ┘
#
#   Key property:  ṽ · w = v × w   for any vector  w ∈ ℝ³.
#
# METHODS
#   tilde(a::VEC3) → MAT3
#   tilde(a::RV3)  → MAT3
#
# USAGE IN SO(3)
#   The skew-symmetric matrix appears in the Rodrigues rotation formula:
#
#     R(ψ) = cos θ · I  +  sin θ · k̃  +  (1 − cos θ) · k⊗k
#
#   and in the linearisation of SO(3) around the identity, where the
#   Lie algebra so(3) is identified with skew-symmetric matrices via `tilde`.
#
# PROPERTIES
#   • tilde(v) is skew-symmetric:  tilde(v)ᵀ = -tilde(v)
#   • tilde(v) has rank 2 (for v ≠ 0) and one zero eigenvalue along v
#   • tilde(v) · v = 0   (the axis is in the null space)
#   • tilde(v) · w = v × w   for all  w
#
# EXTERNAL DEPENDENCIES
#   MAT3, VEC3, RV3    see MAT3.jl, VEC3.jl, RV3.jl
#
# NOTES
#   • The two methods are identical in implementation; separate dispatch
#     avoids the overhead of explicit type conversion between VEC3 and RV3.
# ==============================================================================


"""
    tilde(a::VEC3) → MAT3

Return the 3×3 skew-symmetric (cross-product) matrix of `a`:

    tilde(a) = ┌  0    -a.z   a.y ┐
               │  a.z   0    -a.x │
               └ -a.y   a.x   0   ┘

Such that  `tilde(a) * w = a × w`  for any vector  `w`.
"""
function tilde(a::VEC3)
    return MAT3(
         0.0,  -a.z,   a.y,
         a.z,   0.0,  -a.x,
        -a.y,   a.x,   0.0,
    )
end

"""
    tilde(a::RV3) → MAT3

Return the 3×3 skew-symmetric (cross-product) matrix of the rotation vector `a`.
Identical in computation to the `VEC3` method; provided for dispatch
convenience so that `RV3` arguments need not be explicitly converted.

See `tilde(::VEC3)` for the full definition.
"""
function tilde(a::RV3)
    return MAT3(
         0.0,  -a.z,   a.y,
         a.z,   0.0,  -a.x,
        -a.y,   a.x,   0.0,
    )
end
