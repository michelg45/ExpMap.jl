# ==============================================================================
# DT_SO3.jl — Gâteaux Derivative of the SO(3) Tangent Operator
# ==============================================================================
#
# Computes the 3×3 matrix of the Gâteaux derivative of the map  φ ↦ T(φ)·a
# (or  φ ↦ Tᵀ(φ)·a  when trp = true)  with respect to φ ∈ ℝ³,
# where a ∈ ℝ³ is a fixed vector and T(φ) is the SO(3) tangent operator.
#
# FORMULA (trp = false, default)
#
#   Recall:  T(φ) = I + α(x) φ̂̃ + β(x) (φ̂⊗φ̂ − I)
#
#   where  x = ‖φ‖,  φ̂ = φ/x,  and  α, β  are the scalar coefficients
#   computed by T_functions.  Differentiating  T(φ)·a  w.r.t. φ gives:
#
#   d(T·a)/dφ = [−α·ã + β·DWa] · dφ̂  +  α′·(φ̂×a)⊗φ̂  +  β′·Wa⊗φ̂
#
# FORMULA (trp = true)
#
#   Since  Tᵀ(φ) = I − α(x) φ̃ + β(x) (φ̂⊗φ̂ − I),  flipping the sign
#   of the α term gives:
#
#   d(Tᵀ·a)/dφ = [+α·ã + β·DWa] · dφ̂  −  α′·(φ̂×a)⊗φ̂  +  β′·Wa⊗φ̂
#
#   Setting s = trp ? −1 : +1, both cases unify as:
#
#   (−s·α·ã + β·DWa) · dφ̂  +  s·α′·(φ̂×a)⊗φ̂  +  β′·Wa⊗φ̂
#
# DERIVATION NOTES
#
#   The chain rule is applied to the two scalar-field terms:
#     • α(‖φ‖) · (φ̂×a):  d/dφ = α·d(φ̂×a)/dφ̂·dφ̂ + α′·(φ̂×a)⊗(d‖φ‖/dφ)
#                          with  d(φ̂×a)/dφ̂ = −ã   and  d‖φ‖/dφ = φ̂
#     • β(‖φ‖) · (φ̂⊗φ̂−I)·a:  d/dφ = β·DWa·dφ̂ + β′·Wa⊗φ̂
#   Adding both and grouping the dφ̂ terms gives the formula above.
#
# QUANTITIES
#
#     φ̂   = φ/‖φ‖                        unit rotation axis
#     dφ̂  = (I − φ̂⊗φ̂)/‖φ‖               Gâteaux derivative of (φ → φ̂)  [3×3]
#     φ̂    is also the gradient of ‖φ‖   used as dx in the chain rule   [3×1]
#     Wa   = (φ̂⊗φ̂ − I)·a               β-contribution vector           [3×1]
#     DWa  = φ̂⊗a + (φ̂·a)·I              Gâteaux derivative of (φ̂⊗φ̂)·a  [3×3]
#     ã    = tilde(a)                    skew-symmetric matrix of a     [3×3]
#
# SMALL-ANGLE LIMIT
#
#   For ‖φ‖² < 1e-8:
#     T(φ)  ≈ I − ½φ̃  ⟹  d(T·a)/dφ  ≈  +½·ã   (trp = false)
#     Tᵀ(φ) ≈ I + ½φ̃  ⟹  d(Tᵀ·a)/dφ ≈  −½·ã   (trp = true)
#
# METHOD
#   DT_SO3(φ::RV3, a::VEC3; trp::Bool = false) → MAT3
#
# EXTERNAL DEPENDENCIES
#   T_functions(x; der=1) → (α, α′, β, β′)   scalar coefficients + derivatives
#   tilde(v::VEC3)        → MAT3              skew-symmetric matrix  (ṽ·w = v×w)
#   outerp(u, v)          → MAT3              outer product  u⊗v
#   crossp(u, v)          → VEC3              cross product  u×v
#   dotp(u, v)            → Float64           dot product
#   one(MAT3), VEC3, MAT3
#
# ==============================================================================


"""
    DT_SO3(φ::RV3, a::VEC3; trp::Bool = false) → MAT3

Return the 3×3 matrix of the Gâteaux derivative of `φ ↦ T(φ)·a`
(or of `φ ↦ Tᵀ(φ)·a` when `trp = true`) with respect to `φ`,
where `T(φ)` is the SO(3) tangent operator ([`T_SO3`](@ref)) and `a` is
a fixed vector.

Default (`trp = false`):

    d(T·a)/dφ = [−α·ã + β·DWa] · dφ̂  +  α′·(φ̂×a)⊗φ̂  +  β′·Wa⊗φ̂

Transpose (`trp = true`), using `Tᵀ(φ) = I − α·φ̃ + β·(φ̂⊗φ̂ − I)`:

    d(Tᵀ·a)/dφ = [+α·ã + β·DWa] · dφ̂  −  α′·(φ̂×a)⊗φ̂  +  β′·Wa⊗φ̂

with:
- `φ̂ = φ/‖φ‖` — unit rotation axis
- `dφ̂ = (I − φ̂⊗φ̂)/‖φ‖` — Gâteaux derivative of the normalisation map `φ → φ̂`
- `φ̂` — also the gradient `d‖φ‖/dφ`
- `Wa  = (φ̂⊗φ̂ − I)·a` — β-contribution vector
- `DWa = φ̂⊗a + (φ̂·a)·I` — Gâteaux derivative of `(φ̂⊗φ̂)·a` w.r.t. `φ̂`
- `α, α′, β, β′` — scalar coefficients from [`T_functions`](@ref)

Small-angle limit (`‖φ‖² < 1e-8`): `+½·ã` (`trp = false`) or `-½·ã` (`trp = true`).

This derivative arises in the second-order linearisation of SO(3) kinematics,
in particular in the directional derivative of `T_SO3` needed for tangent
stiffness computations.
"""
function DT_SO3(φ::RV3, a::VEC3; trp::Bool = false)
    x²  = dotp(φ, φ)
    s   = trp ? -1.0 : 1.0
    x² < 1e-8 && return s * 0.5 * tilde(a)   # T ≈ I∓½φ̃  ⟹  d(T·a)/dφ ≈ ±½ã

    I   = one(MAT3)
    x   = sqrt(x²)
    φ̂   = VEC3(φ[1]/x, φ[2]/x, φ[3]/x)   # unit rotation axis

    α, dα, β, dβ = T_functions(x; der=1)   # scalar coefficients and their derivatives

    dφ̂  = (I - outerp(φ̂, φ̂)) / x         # Gâteaux derivative of (φ → φ̂),  [3×3]
    # gradient of ‖φ‖ is φ̂ — substituted directly as the second outerp argument

    Wa  = (outerp(φ̂, φ̂) - I) * a          # β-contribution vector  (φ̂⊗φ̂ − I)·a
    DWa = outerp(φ̂, a) + dotp(φ̂, a) * I   # Gâteaux derivative of (φ̂⊗φ̂)·a w.r.t. φ̂

    return (-s*α * tilde(a) + β * DWa) * dφ̂ +
           s * dα * outerp(crossp(φ̂, a), φ̂) + dβ * outerp(Wa, φ̂)
end
