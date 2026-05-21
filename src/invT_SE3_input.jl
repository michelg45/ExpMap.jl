# ==============================================================================
# invT_SE3_input.jl — Auxiliary quantities for the inverse SE(3) tangent operator
# ==============================================================================
#
# Computes the three quantities needed to assemble invT_SE3(p) = T_SE3(p)⁻¹,
# sharing all intermediate calculations.
#
# INPUT
#   p::VEC6 — 6-vector [u; φ], with translation u = p[1:3] and
#             Cartesian rotation vector φ = p[4:6].
#
# OUTPUT   (MAT3, VEC3, MAT3)
#   invT = T⁻¹(φ)         inverse SO(3) tangent operator
#   Tu   = T(φ)·u          tangent operator applied to the translation part
#   dTu  = DT_SO3(φ, u)    Gâteaux derivative matrix of  φ ↦ T(φ)·u
#
# NOTE
#   Although this function targets the inverse SE(3) operator, the required
#   building blocks are T⁻¹(φ) together with T(φ)·u and DT_SO3(φ,u) — not
#   their T⁻¹ counterparts — because the block formula for T_SE3⁻¹ reads:
#
#     invT_SE3(p) = [ T⁻¹(φ)      T⁻¹(φ)·(−DT_SO3(φ,u)·T⁻¹(φ) + tilde(T(φ)·u)) ]
#                  [ 0₃ₓ₃         T⁻¹(φ)                                           ]
#
#   The shared computation of φ̂, α, β, γ makes the combined evaluation more
#   efficient than calling invT_SO3, T_SO3, and DT_SO3 separately.
#
# SMALL-ANGLE LIMIT   ‖φ‖² < 1e-8
#   invT ≈ I + ½φ̃,   Tu ≈ u − ½(φ×u),   dTu ≈ ½ã
#
# EXTERNAL DEPENDENCIES
#   T_functions(x; der=1), invT_functions(x)
#   tilde, outerp, crossp, dotp, one(MAT3), VEC3, MAT3
#
# ==============================================================================

"""
    invT_SE3_input(p::VEC6) → (invT::MAT3, Tu::VEC3, dTu::MAT3)

Compute the three SO(3) quantities needed to assemble the inverse SE(3)
tangent operator, sharing all intermediate calculations.

Given `p = [u; φ]` with translation `u = p[1:3]` and rotation vector
`φ = p[4:6]`:

- `invT = T⁻¹(φ)`       — inverse SO(3) tangent operator  ([`invT_SO3`](@ref))
- `Tu   = T(φ)·u`        — tangent operator applied to the translation part
- `dTu  = DT_SO3(φ, u)`  — Gâteaux derivative matrix of `φ ↦ T(φ)·u`

The three outputs share the computation of `φ̂`, `α`, `β`, `γ` and their
derivatives, making the combined evaluation more efficient than three
independent calls.
"""
function invT_SE3_input(p::VEC6)
    φ = VEC3(p[4], p[5], p[6])   # rotation vector
    u = VEC3(p[1], p[2], p[3])   # translation vector

    x² = dotp(φ, φ)
    I  = one(MAT3)

    if x² < 1.0e-8
        invT = I + 0.5 * tilde(φ)         # T⁻¹(φ) ≈ I + ½φ̃
        Tu   = u - 0.5 * crossp(φ, u)     # T(φ)·u ≈ u − ½(φ×u)
        dTu  = 0.5 * tilde(u)             # DT_SO3(φ, u) ≈ ½ã
        return invT, Tu, dTu
    end

    x   = sqrt(x²)
    φ̂   = VEC3(φ[1]/x, φ[2]/x, φ[3]/x)        # unit rotation axis

    α, dα, β, dβ = T_functions(x; der=1)
    γ            = invT_functions(x)

    Opp = outerp(φ̂, φ̂)                          # φ̂⊗φ̂
    W   = Opp - I                               # φ̂⊗φ̂ − I
    Wu  = W * u                                 # (φ̂⊗φ̂ − I)·u
    φ̂xu = crossp(φ̂, u)                          # φ̂×u

    invT = I + 0.5 * tilde(φ) + γ * W
    Tu   = u + α * φ̂xu + β * Wu

    dφ̂  = (I - Opp) / x                         # Jacobian of (φ → φ̂),  [3×3]
    DWu = outerp(φ̂, u) + dotp(φ̂, u) * I         # Gâteaux derivative of (φ̂⊗φ̂)·u

    dTu = (-α * tilde(u) + β * DWu) * dφ̂ +
          dα * outerp(φ̂xu, φ̂) + dβ * outerp(Wu, φ̂)

    return invT, Tu, dTu
end
