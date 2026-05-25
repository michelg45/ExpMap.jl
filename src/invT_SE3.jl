# ==============================================================================
# invT_SE3.jl — Inverse SE(3) Tangent Operator
# ==============================================================================
#
# Constructs the 6×6 inverse SE(3) tangent operator T_SE3(p)⁻¹ from either:
#   • a VEC6 input p  (calls invT_SE3_input internally), or
#   • a T_SE3_data struct  (uses precomputed SO(3) data).
#
# For a block upper-triangular operator  T_SE3 = [T  B; 0  T]  the inverse is:
#
#   T_SE3(p)⁻¹ = [ T⁻¹(φ)      T⁻¹(φ)·(−DT_SO3(φ,u)·T⁻¹(φ) + tilde(T(φ)·u)) ]
#                [ 0₃ₓ₃         T⁻¹(φ)                                          ]
#
# which can be written compactly as:
#
#   invT_SE3(p) = [ invT      invT·(−dTu·invT + tilde(Tu)) ]
#                 [ 0₃ₓ₃      invT                         ]
#
# where:
#   invT    = T⁻¹(φ)         inverse SO(3) tangent operator
#   Tu      = T(φ)·u          tangent operator applied to the translation part
#   dTu     = DT_SO3(φ, u)   Gâteaux derivative matrix of  φ ↦ T(φ)·u
#
# sign_p == "−"  VARIANT  (evaluates invT_SE3 at −p)
#
#   Uses the identity T⁻¹(−φ) = (T⁻¹(φ))ᵀ:
#
#   invT_SE3(−p) = [ (T⁻¹)ᵀ      −T⁻¹ · dTu · (T⁻¹)ᵀ ]
#                  [ 0₃ₓ₃          (T⁻¹)ᵀ              ]
#
#   where T⁻¹ = d.invT  and  dTu = d.dTu  (precomputed at +p, reused as-is).
#
# METHODS
#   invT_SE3(p::VEC6)                          → MAT6
#   invT_SE3(d::T_SE3_data, sign_p::String)    → MAT6   ("+" for p, "−" for −p)
#
# EXTERNAL DEPENDENCIES
#   invT_SE3_input, T_SE3_data, tilde, MAT3, MAT6
#
# ==============================================================================

"""
    invT_SE3(p::VEC6) → MAT6
    invT_SE3(d::T_SE3_data, sign_p::String) → MAT6

Construct the 6×6 inverse SE(3) tangent operator `T_SE3(p)⁻¹`.

For `p = [u; φ]` with translation `u = p[1:3]` and rotation vector `φ = p[4:6]`:

    invT_SE3(p) = [ T⁻¹(φ)      T⁻¹(φ)·(−dTu·T⁻¹(φ) + tilde(Tu)) ]
                  [ 0₃ₓ₃         T⁻¹(φ)                             ]

where `invT = T⁻¹(φ)`, `Tu = T(φ)·u`, and `dTu = DT_SO3(φ, u)`.

**First form** — `invT_SE3(p::VEC6)`: calls [`invT_SE3_input`](@ref) internally.

**Second form** — `invT_SE3(d::T_SE3_data, sign_p)`: uses precomputed data
from [`T_SE3_input_data`](@ref).
- `sign_p = "+"` : evaluates `invT_SE3(p)` using `d.invT`, `d.Tu`, `d.dTu`.
- `sign_p = "-"` : evaluates `invT_SE3(−p)` using the identity T⁻¹(−φ) = (T⁻¹(φ))ᵀ:

      invT_SE3(−p) = [ (T⁻¹)ᵀ      −T⁻¹ · dTu · (T⁻¹)ᵀ ]
                     [ 0₃ₓ₃          (T⁻¹)ᵀ              ]

  where `T⁻¹ = d.invT` and `dTu = d.dTu` are reused from the `+p` computation.
"""
function invT_SE3(p::VEC6)
    invT, Tu, dTu = invT_SE3_input(p)
    return MAT6(invT, invT * (-dTu * invT + tilde(Tu)), MAT3(), invT)
end

function invT_SE3(d::T_SE3_data, sign_p::String)
    if sign_p == "+"
        invT = d.invT
        Tu   = d.Tu
        dTu  = d.dTu
        return MAT6(invT, invT * (-dTu * invT + tilde(Tu)), MAT3(), invT)
    else   # sign_p == "-" : evaluate at −p using T⁻¹(−φ) = (T⁻¹(φ))ᵀ
        invT_t = transpose(d.invT)               # (T⁻¹(φ))ᵀ
        return MAT6(invT_t, -(d.invT * d.dTu * invT_t), MAT3(), invT_t)
    end
end
