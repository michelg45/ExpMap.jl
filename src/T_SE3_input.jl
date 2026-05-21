# ==============================================================================
# T_SE3_input.jl — Auxiliary quantities for the SE(3) tangent operator
# ==============================================================================
#
# Computes the three quantities needed to assemble the SE(3) tangent operator,
# sharing all intermediate calculations (φ̂, α, β, α′, β′).
#
# INPUT
#   p::VEC6 — 6-vector [u; φ], with translation u = p[1:3] and
#             Cartesian rotation vector φ = p[4:6].
#
# OUTPUT   (MAT3, VEC3, MAT3)
#   T   = T_SO3(φ)           3×3 SO(3) tangent operator
#   Tu  = T(φ)·u             tangent operator applied to the translation part
#   dTu = DT_SO3(φ, u)       Gâteaux derivative matrix of  φ ↦ T(φ)·u
#
# SMALL-ANGLE LIMIT   ‖φ‖² < 1e-8
#   T   ≈ I − ½φ̃
#   Tu  ≈ u − ½(φ×u)
#   dTu ≈ ½ã
#
# EXTERNAL DEPENDENCIES
#   T_functions(x; der=1), tilde, outerp, crossp, dotp, one(MAT3), VEC3, MAT3
#
# ==============================================================================

"""
    T_SE3_input(p::VEC6) → (T::MAT3, Tu::VEC3, dTu::MAT3)

Compute the three SO(3) quantities needed to assemble the SE(3) tangent
operator, sharing all intermediate calculations.

Given `p = [u; φ]` with translation `u = p[1:3]` and Cartesian rotation
vector `φ = p[4:6]`:

- `T   = T_SO3(φ)`       — 3×3 SO(3) tangent operator  ([`T_SO3`](@ref))
- `Tu  = T(φ)·u`         — tangent operator applied to the translation part
- `dTu = DT_SO3(φ, u)`   — Gâteaux derivative matrix of `φ ↦ T(φ)·u`
  ([`DT_SO3`](@ref))

The combined evaluation is more efficient than three independent calls because
`φ̂`, `α`, `β`, `α′`, `β′` and the outer product `φ̂⊗φ̂` are computed once
and shared across all three outputs.
"""
function T_SE3_input(p::VEC6)
    φ = VEC3(p[4], p[5], p[6])   # rotation vector
    u = VEC3(p[1], p[2], p[3])   # translation vector

    x² = dotp(φ, φ)
    I  = one(MAT3)

    if x² < 1.0e-8
        T   = I - 0.5 * tilde(φ)         # T(φ) ≈ I − ½φ̃
        Tu  = u - 0.5 * crossp(φ, u)     # T(φ)·u ≈ u − ½(φ×u)
        dTu = 0.5 * tilde(u)             # DT_SO3(φ, u) ≈ ½ã
        return T, Tu, dTu
    end

    x   = sqrt(x²)
    φ̂   = VEC3(φ[1]/x, φ[2]/x, φ[3]/x)

    α, dα, β, dβ = T_functions(x; der=1)

    Opp = outerp(φ̂, φ̂)               # φ̂⊗φ̂                              [3×3]
    W   = Opp - I                     # φ̂⊗φ̂ − I                          [3×3]
    Wa  = W * u                       # (φ̂⊗φ̂ − I)·u                     [3×1]
    φ̂xu = crossp(φ̂, u)               # φ̂ × u                            [3×1]

    T   = I + α * tilde(φ̂) + β * W
    Tu  = u + α * φ̂xu + β * Wa
    dTu = (-α * tilde(u) + β * (outerp(φ̂, u) + dotp(φ̂, u) * I)) * ((I - Opp) / x) +
          dα * outerp(φ̂xu, φ̂) + dβ * outerp(Wa, φ̂)

    return T, Tu, dTu
end
