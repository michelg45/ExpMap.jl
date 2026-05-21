# ==============================================================================
# Adj_SE3.jl — Adjoint Representation of SE(3) and se(3) Lie Bracket
# ==============================================================================
#
# Provides four overloads implementing:
#   1. the group adjoint  Ad(H)  of SE(3),
#   2. the matrix-free adjoint action  Ad(H)·a,
#   3. the Lie-algebra small adjoint  ad(a) (= sk_SE3(a)),
#   4. the Lie bracket  [a, b]  computed directly on vectors.
#
# ── GROUP ADJOINT  Ad(H) ────────────────────────────────────────────────────
#
#   For H = (x, φ) ∈ SE(3) with rotation matrix R = R_SO3(φ):
#
#     Ad(H) = [R    x̃·R]   ∈ ℝ⁶ˣ⁶
#             [0    R  ]
#
#   where  x̃ = tilde(x).  For any twist a = (a_u; a_ω):
#
#     (Ad(H)·a)[1:3] = R·a_u  +  x × (R·a_ω)
#     (Ad(H)·a)[4:6] = R·a_ω
#
#   Geometrically, Ad(H) changes the reference frame of a twist: if a is
#   expressed in body frame, Ad(H)·a is the same twist in spatial frame.
#
# ── SMALL ADJOINT  ad(a) ────────────────────────────────────────────────────
#
#   For a = (u; ω) ∈ se(3):
#
#     ad(a) = [ω̃   ũ]   ∈ ℝ⁶ˣ⁶
#             [0    ω̃]
#
#   ad(a) equals sk_SE3(a) and satisfies  ad(a)·b = [a, b].
#
# ── LIE BRACKET  [a, b] ─────────────────────────────────────────────────────
#
#   For a = (u_a; ω_a) and b = (u_b; ω_b):
#
#     [a, b][1:3] = ω_a × u_b  +  u_a × ω_b
#     [a, b][4:6] = ω_a × ω_b
#
# METHODS
#   Adj_SE3(H::NodeFrame)           → MAT6   full 6×6 adjoint matrix Ad(H)
#   Adj_SE3(H::NodeFrame, a::VEC6) → VEC6   adjoint action Ad(H)·a (no matrix)
#   Adj_SE3(a::VEC6)                → MAT6   small adjoint ad(a)
#   Adj_SE3(a::VEC6, b::VEC6)      → VEC6   Lie bracket [a, b] (no matrix)
#
# EXTERNAL DEPENDENCIES
#   R_SO3(psi::RV3)           → MAT3   rotation matrix R(φ)
#   R_SO3(psi::RV3, a::VEC3) → VEC3   matrix-free rotation R(φ)·a
#   tilde(v::VEC3)            → MAT3   skew-symmetric matrix
#   crossp(u::VEC3, v::VEC3) → VEC3   cross product  u × v
#   VEC3, VEC6, MAT3, MAT6, NodeFrame
#
# ==============================================================================


"""
    Adj_SE3(H::NodeFrame) → MAT6

Return the 6×6 adjoint matrix `Ad(H)` for the SE(3) element `H = (x, φ)`:

    Ad(H) = [R    x̃·R]
            [0    R  ]

where `R = R_SO3(H.phi)` and `x̃ = tilde(H.x)`.

`Ad(H)` transforms a twist from body frame to spatial frame: for a twist
`a = (a_u; a_ω)`, `Ad(H)·a` gives the spatially-expressed equivalent twist.
"""
function Adj_SE3(H::NodeFrame)
    R = R_SO3(H.phi)                              # rotation matrix R = R(φ)
    return MAT6(R, tilde(H.x) * R, zero(MAT3), R) # [R  x̃·R ; 0  R]
end


"""
    Adj_SE3(H::NodeFrame, a::VEC6) → VEC6

Apply the adjoint `Ad(H)` to the twist `a = (a_u; a_ω)` **without forming
the full 6×6 matrix**:

    (Ad(H)·a)[1:3] = R·a_u  +  x × (R·a_ω)
    (Ad(H)·a)[4:6] = R·a_ω

where `R = R_SO3(H.phi)` and `x = H.x`.

!!! tip "Efficiency"
    Prefer `Adj_SE3(H, a)` over `Adj_SE3(H) * a` when only the image
    vector is needed — it avoids allocating the full 6×6 matrix.
"""
function Adj_SE3(H::NodeFrame, a::VEC6)
    a_u     = VEC3(a[1], a[2], a[3])           # translational part of a
    a_omega = VEC3(a[4], a[5], a[6])           # angular part of a
    Ra_omega = R_SO3(H.phi, a_omega)           # R·a_ω  (matrix-free)
    Ra_u     = R_SO3(H.phi, a_u)              # R·a_u  (matrix-free)
    return VEC6(Ra_u + crossp(H.x, Ra_omega), Ra_omega)  # [R·a_u + x×(R·a_ω); R·a_ω]
end


"""
    Adj_SE3(a::VEC6) → MAT6

Return the 6×6 small-adjoint matrix `ad(a)` for the se(3) element
`a = (u; ω)`:

    ad(a) = [ω̃   ũ]
            [0    ω̃]

`ad(a)·b = [a, b]` is the Lie bracket of `a` and `b` in se(3).
This matrix is identical to [`sk_SE3`](@ref) applied to `a`.
"""
function Adj_SE3(a::VEC6)
    a_u     = VEC3(a[1], a[2], a[3])   # translational part  u
    a_omega = VEC3(a[4], a[5], a[6])   # angular part  ω
    return MAT6(tilde(a_omega), tilde(a_u), zero(MAT3), tilde(a_omega))  # [ω̃  ũ; 0  ω̃]
end


"""
    Adj_SE3(a::VEC6, b::VEC6) → VEC6

Return the Lie bracket `[a, b]` of two se(3) elements
`a = (u_a; ω_a)` and `b = (u_b; ω_b)`, **without forming the 6×6 matrix**:

    [a, b][1:3] = ω_a × u_b  +  u_a × ω_b
    [a, b][4:6] = ω_a × ω_b

Equivalent to `Adj_SE3(a) * b` but avoids allocating the 6×6 matrix.
"""
function Adj_SE3(a::VEC6, b::VEC6)
    a_u     = VEC3(a[1], a[2], a[3])
    a_omega = VEC3(a[4], a[5], a[6])
    b_u     = VEC3(b[1], b[2], b[3])
    b_omega = VEC3(b[4], b[5], b[6])
    return VEC6(crossp(a_omega, b_u) + crossp(a_u, b_omega),   # ω_a×u_b + u_a×ω_b
                crossp(a_omega, b_omega))                        # ω_a×ω_b
end
