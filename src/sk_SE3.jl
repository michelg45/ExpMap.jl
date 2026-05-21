# ==============================================================================
# sk_SE3.jl — Hat Map (Skew Operator) for SE(3)
# ==============================================================================
#
# Implements the hat map for SE(3), mapping a 6-component twist vector
# v = (u; ω) ∈ ℝ⁶ to its 6×6 block-skew matrix representation in se(3).
#
# The twist vector v is partitioned as:
#
#   v[1:3] = u ∈ ℝ³   translational part
#   v[4:6] = ω ∈ ℝ³   angular (rotational) part
#
# HAT MAP FORMULA
#
#   sk_SE3(v) = [ω̃   ũ]    ∈ ℝ⁶ˣ⁶
#               [0    ω̃]
#
# where  ω̃ = tilde(ω)  and  ũ = tilde(u)  are 3×3 skew-symmetric matrices.
#
# ACTION ON A TWIST
#
#   For any a = (a_u; a_ω) ∈ ℝ⁶:
#
#   sk_SE3(v) · a = [ω × a_u + u × a_ω]  =  [v, a]   (Lie bracket of se(3))
#                   [ω × a_ω           ]
#
# RELATION TO Adj_SE3
#
#   sk_SE3(v)  ≡  Adj_SE3(v::VEC6)   (same matrix, different semantics:
#   sk_SE3 emphasises the hat map; Adj_SE3 emphasises the adjoint action)
#
# METHOD
#   sk_SE3(v::VEC6) → MAT6
#
# EXTERNAL DEPENDENCIES
#   tilde(v::VEC3)   → MAT3   skew-symmetric (cross-product) matrix
#   VEC3, VEC6, MAT3, MAT6   see VEC3.jl, VEC6.jl, MAT3.jl, MAT6.jl
#
# ==============================================================================


"""
    sk_SE3(v::VEC6) → MAT6

Return the 6×6 block-skew matrix that represents the SE(3) Lie algebra
element associated with the twist `v = (u; ω)`:

    sk_SE3(v) = [ω̃   ũ]
                [0    ω̃]

where `ω = v[4:6]` (angular part) and `u = v[1:3]` (translational part).

For any twist `a = (a_u; a_ω)` the matrix action recovers the Lie bracket:

    sk_SE3(v) · a = [ω × a_u + u × a_ω]  =  [v, a]
                    [ω × a_ω           ]
"""
function sk_SE3(v::VEC6)
    omega = VEC3(v[4], v[5], v[6])   # angular part  ω = v[4:6]
    u     = VEC3(v[1], v[2], v[3])   # translational part  u = v[1:3]
    return MAT6(tilde(omega), tilde(u), zero(MAT3), tilde(omega))  # [ω̃  ũ; 0  ω̃]
end
