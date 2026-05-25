# ==============================================================================
# T_SE3.jl — SE(3) Tangent Operator
# ==============================================================================
#
# Constructs the 6×6 SE(3) tangent operator from either:
#   • a VEC6 input p  (calls T_SE3_input internally), or
#   • a T_SE3_data struct  (uses precomputed SO(3) data).
#
# For parameter vector p = [u; φ] with translation u ∈ ℝ³ and rotation
# vector φ ∈ ℝ³, the SE(3) tangent operator reads:
#
#   T_SE3(p) = [ T(φ)      dTu − tilde(Tu) · T(φ) ]
#              [ 0₃ₓ₃      T(φ)                    ]
#
# where:
#   T(φ)    = T_SO3(φ)           SO(3) tangent operator
#   Tu      = T(φ)·u             translation contribution
#   dTu     = DT_SO3(φ, u)       Gâteaux derivative matrix of φ ↦ T(φ)·u
#
# ACTION ON A VARIATION  [δu; δφ]:
#
#   T_SE3(p) · [δu; δφ] = [ T·δu + (dTu − tilde(Tu)·T)·δφ ]
#                          [ T·δφ                            ]
#
# sign_p == "−"  VARIANT  (evaluates T_SE3 at −p = [−u; −φ])
#
#   Uses the identity T(−φ) = T(φ)ᵀ:
#     T        ← T(φ)ᵀ
#     Tu       ← −T(φ)ᵀ·u
#     dTu[:,i] ← (∂T/∂φᵢ)ᵀ·u
#   then applies the same block formula.
#
# METHODS
#   T_SE3(p::VEC6)                          → MAT6   (forward, VEC6 input)
#   T_SE3(d::T_SE3_data, sign_p::String)    → MAT6   ("+" for p, "−" for −p)
#
# EXTERNAL DEPENDENCIES
#   T_SE3_input, T_SE3_data, tilde, MAT3, MAT6, VEC3
#
# ==============================================================================

"""
    T_SE3(p::VEC6) → MAT6
    T_SE3(d::T_SE3_data, sign_p::String) → MAT6

Construct the 6×6 SE(3) tangent operator.

For `p = [u; φ]` with translation `u = p[1:3]` and rotation vector `φ = p[4:6]`:

    T_SE3(p) = [ T(φ)      dTu − tilde(Tu) · T(φ) ]
               [ 0₃ₓ₃      T(φ)                    ]

where `T(φ) = T_SO3(φ)`, `Tu = T(φ)·u`, and `dTu = DT_SO3(φ, u)`.

**First form** — `T_SE3(p::VEC6)`: calls [`T_SE3_input`](@ref) to obtain
`T`, `Tu`, `dTu` before assembling the block matrix.

**Second form** — `T_SE3(d::T_SE3_data, sign_p)`: uses precomputed data
from [`T_SE3_input_data`](@ref).
- `sign_p = "+"` : evaluates `T_SE3(p)` using `d.T`, `d.Tu`, `d.dTu`.
- `sign_p = "-"` : evaluates `T_SE3(−p)` using the identity T(−φ) = T(φ)ᵀ,
  substituting `T ← T(φ)ᵀ`, `Tu ← −T(φ)ᵀ·u`, `dTu[:,i] ← (∂T/∂φᵢ)ᵀ·u`.
"""
function T_SE3(p::VEC6)
    T, Tu, dTu = T_SE3_input(p)
    return MAT6(T, dTu - tilde(Tu) * T, MAT3(), T)
end

function T_SE3(d::T_SE3_data, sign_p::String)
    if sign_p == "+"
        T   = d.T
        Tu  = d.Tu
        dTu = d.dTu
    else   # sign_p == "-" : evaluate at −p using T(−φ) = T(φ)ᵀ
        u   = VEC3(d.p[1], d.p[2], d.p[3])
        T   = transpose(d.T)
        Tu  = -(T * u)                           # T(−φ)·(−u) = −T(φ)ᵀ·u
        dTu = MAT3(transpose(d.dT[1]) * u,
                   transpose(d.dT[2]) * u,
                   transpose(d.dT[3]) * u)       # (∂T/∂φᵢ)ᵀ·u
    end
    return MAT6(T, dTu - tilde(Tu) * T, MAT3(), T)
end
