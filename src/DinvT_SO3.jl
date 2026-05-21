# ==============================================================================
# DinvT_SO3.jl — Gâteaux Derivative of the Inverse SO(3) Tangent Operator
# ==============================================================================
#
# Computes the 3×3 matrix of the Gâteaux derivative of the map  φ ↦ T⁻¹(φ)·a
# (or  φ ↦ (T⁻¹)ᵀ(φ)·a  when trp = true)  with respect to φ ∈ ℝ³,
# where a ∈ ℝ³ is a fixed vector and T⁻¹(φ) is the inverse tangent operator.
#
# FORMULA (trp = false, default)
#
#   Recall:  T⁻¹(φ) = I + ½ φ̃ + γ(x) (φ̂⊗φ̂ − I)
#
#   where  x = ‖φ‖,  φ̂ = φ/x,  and  γ(x) = 1 − (x/2) cot(x/2)  is the
#   scalar coefficient computed by invT_functions.
#
#   The key structural difference from DT_SO3 is that the linear term uses
#   φ (not φ̂), so its Gâteaux derivative is the constant matrix −½ã:
#
#     d(½ φ̃·a)/dφ = d(½ φ×a)/dφ = −½ ã
#
#   The γ-term is differentiated exactly as in DT_SO3:
#
#     d(γ·(φ̂⊗φ̂−I)·a)/dφ = γ·DWa·dφ̂ + γ′·Wa⊗φ̂
#
#   Combining:
#
#   d(T⁻¹·a)/dφ = −½·ã + γ·DWa·dφ̂ + γ′·Wa⊗φ̂
#
# FORMULA (trp = true)
#
#   Since  (T⁻¹)ᵀ(φ) = I − ½ φ̃ + γ(x) (φ̂⊗φ̂ − I),
#   the sign of the constant term flips:
#
#   d((T⁻¹)ᵀ·a)/dφ = +½·ã + γ·DWa·dφ̂ + γ′·Wa⊗φ̂
#
#   Setting s = trp ? −1 : +1, both cases unify as:
#
#   d(T̂⁻¹·a)/dφ = −s·½·ã + γ·DWa·dφ̂ + γ′·Wa⊗φ̂
#
# DERIVATION NOTES
#
#   ½ φ̃·a = ½ (φ × a)  is linear in φ (unnormalized), so:
#     d(½ φ×a)/dφ · δφ = ½ (δφ×a) = −½ (a×δφ) = −½ ã·δφ   ⟹   matrix = −½ã
#   No dφ̂ contribution arises — contrast with α(x)·(φ̂×a) in DT_SO3.
#
#   The γ-term follows the standard chain rule:
#     d(γ·(φ̂⊗φ̂)·a)/dφ = γ·DWa·dφ̂ + γ′·Wa⊗φ̂
#   where  DWa = ∂((φ̂⊗φ̂)·a)/∂φ̂  and  dφ̂ = (I−φ̂⊗φ̂)/x.
#
# QUANTITIES
#
#     φ̂   = φ/‖φ‖                        unit rotation axis
#     dφ̂  = (I − φ̂⊗φ̂)/‖φ‖               Gâteaux derivative of (φ → φ̂)  [3×3]
#     φ̂    is also the gradient of ‖φ‖   used as dx in the chain rule   [3×1]
#     Wa   = (φ̂⊗φ̂ − I)·a               γ-contribution vector           [3×1]
#     DWa  = φ̂⊗a + (φ̂·a)·I              Gâteaux derivative of (φ̂⊗φ̂)·a  [3×3]
#     ã    = tilde(a)                    skew-symmetric matrix of a     [3×3]
#
# SMALL-ANGLE LIMIT
#
#   For ‖φ‖² < 1e-8:
#     T⁻¹(φ)  ≈ I + ½φ̃   ⟹   d(T⁻¹·a)/dφ  ≈  −½·ã   (trp = false)
#     (T⁻¹)ᵀ  ≈ I − ½φ̃   ⟹   d((T⁻¹)ᵀ·a)/dφ ≈ +½·ã   (trp = true)
#
# METHOD
#   DinvT_SO3(φ::RV3, a::VEC3; trp::Bool = false) → MAT3
#
# EXTERNAL DEPENDENCIES
#   invT_functions(x; der=1) → (γ, γ′)   scalar coefficient + derivative
#   tilde(v::VEC3)           → MAT3       skew-symmetric matrix  (ṽ·w = v×w)
#   outerp(u, v)             → MAT3       outer product  u⊗v
#   dotp(u, v)               → Float64    dot product
#   one(MAT3), VEC3, MAT3
#
# ==============================================================================


"""
    DinvT_SO3(φ::RV3, a::VEC3; trp::Bool = false) → MAT3

Return the 3×3 matrix of the Gâteaux derivative of `φ ↦ T⁻¹(φ)·a`
(or of `φ ↦ (T⁻¹)ᵀ(φ)·a` when `trp = true`) with respect to `φ`,
where `T⁻¹(φ)` is the inverse SO(3) tangent operator ([`invT_SO3`](@ref))
and `a` is a fixed vector.

Default (`trp = false`):

    d(T⁻¹·a)/dφ = −½·ã + γ·DWa·dφ̂ + γ′·Wa⊗φ̂

Transpose (`trp = true`), using `(T⁻¹)ᵀ(φ) = I − ½φ̃ + γ(φ̂⊗φ̂ − I)`:

    d((T⁻¹)ᵀ·a)/dφ = +½·ã + γ·DWa·dφ̂ + γ′·Wa⊗φ̂

with:
- `φ̂ = φ/‖φ‖` — unit rotation axis
- `dφ̂ = (I − φ̂⊗φ̂)/‖φ‖` — Gâteaux derivative of the normalisation map `φ → φ̂`
- `φ̂` — also the gradient `d‖φ‖/dφ`
- `Wa  = (φ̂⊗φ̂ − I)·a` — γ-contribution vector
- `DWa = φ̂⊗a + (φ̂·a)·I` — Gâteaux derivative of `(φ̂⊗φ̂)·a` w.r.t. `φ̂`
- `γ, γ′` — scalar coefficient from [`invT_functions`](@ref)

Note: the `−½·ã` term is a **constant matrix** (independent of `φ̂`), because
`½φ̃·a = ½(φ×a)` is linear in `φ` (not in `φ̂`).

Small-angle limit (`‖φ‖² < 1e-8`): `−½·ã` (`trp = false`) or `+½·ã` (`trp = true`).

This derivative arises in the second-order linearisation of SO(3) kinematics,
for instance in the linearisation of strain measures expressed via `T⁻¹`.
"""
function DinvT_SO3(φ::RV3, a::VEC3; trp::Bool = false)
    x²  = dotp(φ, φ)
    s   = trp ? -1.0 : 1.0
    x² < 1e-8 && return -s * 0.5 * tilde(a)   # T⁻¹ ≈ I + ½φ̃  ⟹  d(T⁻¹·a)/dφ ≈ −½ã

    I   = one(MAT3)
    x   = sqrt(x²)
    φ̂   = VEC3(φ[1]/x, φ[2]/x, φ[3]/x)   # unit rotation axis

    γ, dγ = invT_functions(x; der=1)       # scalar coefficient and its derivative

    dφ̂  = (I - outerp(φ̂, φ̂)) / x         # Gâteaux derivative of (φ → φ̂),  [3×3]
    # gradient of ‖φ‖ is φ̂ — substituted directly as the second outerp argument

    Wa  = (outerp(φ̂, φ̂) - I) * a          # γ-contribution vector  (φ̂⊗φ̂ − I)·a
    DWa = outerp(φ̂, a) + dotp(φ̂, a) * I   # Gâteaux derivative of (φ̂⊗φ̂)·a w.r.t. φ̂

    return -s * 0.5 * tilde(a) + γ * DWa * dφ̂ + dγ * outerp(Wa, φ̂)
end
