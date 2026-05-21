# ==============================================================================
# log_SE3.jl — Logarithmic Map of SE(3): NodeFrame → Motion Parameters
# ==============================================================================
#
# Implements the logarithmic map of the Special Euclidean group SE(3), mapping
# a NodeFrame H = (x, φ) ∈ SE(3) to the motion-parameter vector
# p = (u; φ) ∈ ℝ⁶.
#
# LOGARITHMIC MAP FORMULA
#
#   p[4:6] = φ               (CRV, unchanged)
#   p[1:3] = T⁻¹(-φ) · x   (translational parameter)
#
# where T⁻¹(-φ) is the inverse of the SO(3) tangent operator at −φ:
#
#   T⁻¹(φ) = I + ½ φ̃ + γ(‖φ‖) (φ̂⊗φ̂ − I)
#
# with  γ(x) = 1 − (x/2) cot(x/2).
#
# SMALL-ANGLE LIMIT  (handled transparently by invT_SO3)
#
#   For ‖φ‖ < 1e-4:  T⁻¹(-φ) ≈ I + ½φ̃,  so  u ≈ x + ½ φ × x
#   For ‖φ‖ = 0  :  T⁻¹(0)  = I,          so  u = x  (pure translation)
#
# ROUND-TRIP IDENTITY
#
#   log_SE3(exp_SE3(p)) = p   for all p ∈ ℝ⁶
#   exp_SE3(log_SE3(H)) = H   for all H ∈ SE(3)
#
# METHOD
#   log_SE3(H::NodeFrame) → VEC6
#
# EXTERNAL DEPENDENCIES
#   invT_SO3(phi::RV3, a::VEC3) → VEC3   inverse tangent operator  T⁻¹(φ)·a
#   VEC3, RV3, VEC6, NodeFrame            see VEC3.jl, RV3.jl, VEC6.jl, NodeFrame.jl
#
# ==============================================================================


"""
    log_SE3(H::NodeFrame) → VEC6

Return the motion-parameter vector `p = (u; φ)` corresponding to the SE(3)
element `H`, via the logarithmic map:

    p[4:6] = H.phi                (CRV encoding orientation)
    p[1:3] = T⁻¹(-H.phi) · H.x  (translational parameter)

The inverse tangent operator `T⁻¹(-φ)` recovers the translational parameter
`u` from the position vector `x`.  The map is the exact inverse of
[`exp_SE3`](@ref): if `H = exp_SE3(p)` then `log_SE3(H) = p`.

For `‖φ‖ = 0` (no rotation), `T⁻¹(0) = I` and `p[1:3] = H.x` exactly.

Small-angle handling (`‖φ‖ < 1e-4`) is delegated to [`invT_SO3`](@ref).
"""
function log_SE3(H::NodeFrame)
    return VEC6(invT_SO3(-H.phi, H.x), H.phi)  # u = T⁻¹(-φ)·x,  φ = H.phi
end
