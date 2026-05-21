# ==============================================================================
# exp_SE3.jl — Exponential Map of SE(3): Motion Parameters → NodeFrame
# ==============================================================================
#
# Implements the exponential map of the Special Euclidean group SE(3), mapping
# a 6-component motion-parameter vector p = (u; φ) ∈ ℝ⁶ to a NodeFrame
# H = (x, φ) ∈ SE(3).
#
# The motion-parameter vector p is partitioned as:
#
#   p[1:3] = u ∈ ℝ³   translational parameter  (reduced translation)
#   p[4:6] = φ ∈ ℝ³   Cartesian Rotation Vector (CRV) encoding orientation
#
# EXPONENTIAL MAP FORMULA
#
#   H.phi = φ              (rotation CRV, unchanged)
#   H.x   = T(-φ) · u     (position vector)
#
# where T(-φ) is the SO(3) tangent operator evaluated at −φ:
#
#   T(φ) = I + α(‖φ‖) φ̃/‖φ‖ + β(‖φ‖) (φ̂⊗φ̂ − I)
#
# with  α(x) = (cos x − 1)/x  and  β(x) = 1 − sin(x)/x.
#
# SMALL-ANGLE LIMIT  (handled transparently by T_SO3)
#
#   For ‖φ‖ < 1e-4:  T(-φ) ≈ I − ½φ̃,  so  x ≈ u − ½ φ × u
#   For ‖φ‖ = 0  :  T(0)  = I,          so  x = u  (pure translation)
#
# ROUND-TRIP IDENTITY
#
#   exp_SE3(log_SE3(H)) = H   for all H ∈ SE(3)
#   log_SE3(exp_SE3(p)) = p   for all p ∈ ℝ⁶
#
# METHOD
#   exp_SE3(p::VEC6) → NodeFrame
#
# EXTERNAL DEPENDENCIES
#   T_SO3(phi::RV3, a::VEC3) → VEC3   tangent operator action  T(φ)·a
#   VEC3, RV3, VEC6, NodeFrame         see VEC3.jl, RV3.jl, VEC6.jl, NodeFrame.jl
#
# ==============================================================================


"""
    exp_SE3(p::VEC6) → NodeFrame

Return the SE(3) group element `H` corresponding to the motion-parameter
vector `p = (u; φ)`, via the exponential map:

    H.phi = φ            (CRV encoding orientation)
    H.x   = T(-φ) · u   (position vector)

The vector `p` is partitioned as:
- `p[1:3]` → translational parameter `u ∈ ℝ³`
- `p[4:6]` → Cartesian Rotation Vector `φ ∈ ℝ³`

The position is obtained by applying the SO(3) tangent operator `T(-φ)` to
`u`.  For `‖φ‖ = 0` (no rotation), `T(0) = I` and `H.x = u` exactly.

Small-angle handling (`‖φ‖ < 1e-4`) is delegated to [`T_SO3`](@ref).
"""
function exp_SE3(p::VEC6)
    phi = RV3(p[4], p[5], p[6])               # rotation CRV φ = p[4:6]
    x   = T_SO3(-phi, VEC3(p[1], p[2], p[3])) # position: x = T(-φ)·u
    return NodeFrame(x, phi)
end
