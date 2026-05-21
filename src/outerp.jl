# ==============================================================================
# outerp.jl — Outer (Dyadic) Product of Two 3D Vectors
# ==============================================================================
#
# Provides two methods of `outerp` computing the outer product (also called
# the dyadic or tensor product) of two 3D vectors, returning a 3×3 matrix.
#
# DEFINITION
#   Given vectors  u, v ∈ ℝ³,  the outer product is the rank-1 matrix:
#
#       (u ⊗ v)ᵢⱼ = uᵢ · vⱼ
#
#   ┌ u.x·v.x   u.x·v.y   u.x·v.z ┐
#   │ u.y·v.x   u.y·v.y   u.y·v.z │
#   └ u.z·v.x   u.z·v.y   u.z·v.z ┘
#
# METHODS
#   outerp(v1::VEC3, v2::VEC3) → MAT3
#   outerp(v1::RV3,  v2::RV3)  → MAT3
#
# USAGE IN SO(3)
#   The outer product appears in the Rodrigues rotation formula:
#
#     R(ψ) = cos θ · I  +  sin θ · k̃  +  (1 − cos θ) · outerp(k, k)
#
#   where  k = ψ/‖ψ‖  is the unit rotation axis  (see R_SO3.jl).
#
# EXTERNAL DEPENDENCIES
#   MAT3, VEC3, RV3    see MAT3.jl, VEC3.jl, RV3.jl
#
# NOTES
#   • outerp(u, v) ≠ outerp(v, u) in general (the result is only symmetric
#     when u = v, giving a rank-1 projection matrix).
#   • The two methods are identical in implementation; separate dispatch
#     allows both VEC3 and RV3 arguments without explicit conversion.
# ==============================================================================


"""
    outerp(v1::VEC3, v2::VEC3) → MAT3

Return the outer (dyadic) product  `v1 ⊗ v2`,  the 3×3 rank-1 matrix
whose `(i,j)` entry is `v1[i] · v2[j]`:

    ┌ v1.x·v2.x   v1.x·v2.y   v1.x·v2.z ┐
    │ v1.y·v2.x   v1.y·v2.y   v1.y·v2.z │
    └ v1.z·v2.x   v1.z·v2.y   v1.z·v2.z ┘
"""
function outerp(v1::VEC3, v2::VEC3)
    return MAT3(
        v1.x * v2.x,  v1.x * v2.y,  v1.x * v2.z,
        v1.y * v2.x,  v1.y * v2.y,  v1.y * v2.z,
        v1.z * v2.x,  v1.z * v2.y,  v1.z * v2.z,
    )
end

"""
    outerp(v1::RV3, v2::RV3) → MAT3

Return the outer (dyadic) product  `v1 ⊗ v2`  for two rotation vectors.
Identical in computation to the `VEC3` method; provided for dispatch
convenience so that `RV3` arguments need not be explicitly converted.

See `outerp(::VEC3, ::VEC3)` for the full definition.
"""
function outerp(v1::RV3, v2::RV3)
    return MAT3(
        v1.x * v2.x,  v1.x * v2.y,  v1.x * v2.z,
        v1.y * v2.x,  v1.y * v2.y,  v1.y * v2.z,
        v1.z * v2.x,  v1.z * v2.y,  v1.z * v2.z,
    )
end
