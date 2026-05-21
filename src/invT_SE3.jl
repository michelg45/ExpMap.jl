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
# TRANSPOSE VARIANT  (trp = true, requires T_SE3_data input)
#
#   Substitutes  invT → (T⁻¹)ᵀ,  Tu → −(T⁻¹)ᵀ·u,  dTu[:,i] → (∂T/∂φᵢ)ᵀ·u
#   into the same block formula.
#
# METHODS
#   invT_SE3(p::VEC6)                         → MAT6   (forward, VEC6 input)
#   invT_SE3(d::T_SE3_data; trp::Bool=false)  → MAT6   (from T_SE3_data)
#
# EXTERNAL DEPENDENCIES
#   invT_SE3_input, T_SE3_data, tilde, MAT3, MAT6, VEC3
#
# ==============================================================================

"""
    invT_SE3(p::VEC6) → MAT6
    invT_SE3(d::T_SE3_data; trp::Bool = false) → MAT6

Construct the 6×6 inverse SE(3) tangent operator `T_SE3(p)⁻¹`.

For `p = [u; φ]` with translation `u = p[1:3]` and rotation vector `φ = p[4:6]`:

    invT_SE3(p) = [ T⁻¹(φ)      T⁻¹(φ)·(−dTu·T⁻¹(φ) + tilde(Tu)) ]
                  [ 0₃ₓ₃         T⁻¹(φ)                             ]

where `invT = T⁻¹(φ)`, `Tu = T(φ)·u`, and `dTu = DT_SO3(φ, u)`.

**First form** — `invT_SE3(p::VEC6)`: calls [`invT_SE3_input`](@ref) internally.

**Second form** — `invT_SE3(d::T_SE3_data; trp = false)`: uses precomputed data
from [`T_SE3_input_data`](@ref).  With `trp = true`, replaces:
- `invT`       →  `(T⁻¹)ᵀ(φ)`
- `Tu`          →  `−(T⁻¹)ᵀ(φ)·u`
- `dTu[:,i]`    →  `(∂T/∂φᵢ)ᵀ·u`

and applies the same block formula to the modified triple.
"""
function invT_SE3(p::VEC6)
    invT, Tu, dTu = invT_SE3_input(p)
    return MAT6(invT, invT * (-dTu * invT + tilde(Tu)), MAT3(), invT)
end

function invT_SE3(d::T_SE3_data; trp::Bool = false)
    if !trp
        invT = d.invT
        Tu   = d.Tu
        dTu  = d.dTu
    else
        u    = VEC3(d.p[1], d.p[2], d.p[3])
        invT = transpose(d.invT)
        Tu   = -(invT * u)                       # Tu = −(T⁻¹)ᵀ(φ)·u
        dTu  = MAT3(transpose(d.dT[1]) * u,
                    transpose(d.dT[2]) * u,
                    transpose(d.dT[3]) * u)      # dTu[:,i] = (∂T/∂φᵢ)ᵀ·u
    end
    return MAT6(invT, invT * (-dTu * invT + tilde(Tu)), MAT3(), invT)
end
