# ==============================================================================
# DinvT_SE3.jl — Gâteaux Derivative of the Inverse SE(3) Tangent Operator
# ==============================================================================
#
# Computes the 6×6 matrix of the Gâteaux derivative of  p ↦ invT_SE3(p)·f
# (or  p ↦ invT_SE3(p)ᵀ·f)  with respect to the full parameter vector p = [u; φ].
#
# Recall the block structure of the inverse operator:
#
#   invT_SE3(p) = [ invT    B  ]     B = invT·(−dTu·invT + tilde(Tu))
#                 [  0    invT ]
#
# where  invT = T⁻¹(φ),  Tu = T(φ)·u,  dTu = DT_SO3(φ,u).
#
# INPUT
#   a       :: T_SE3_data   precomputed data from T_SE3_input_data
#   f       :: VEC6         fixed 6-vector; split as f_u = f[1:3], f_φ = f[4:6]
#   trp     :: Bool = false transpose flag (true → differentiate invT_SE3ᵀ·f)
#   sign_p  :: String = "+" operator sign ("+" → forward invT_SE3; "-" → uses invTᵀ)
#
# OUTPUT   MAT6   — the 6×6 Gâteaux derivative matrix written in block form:
#
#   [ D_uu  D_uφ ]     rows 1–3: derivative of upper-block output
#   [ D_φu  D_φφ ]     rows 4–6: derivative of lower-block output
#
# FOUR CASES   (trp) × (sign_p)
#
#   (!trp, "+")  D_p(invT_SE3(p)·f)              forward, no transpose
#   ( trp, "+")  D_p(invT_SE3(p)ᵀ·f)            forward, transposed
#   (!trp, "-")  D_p(invT_SE3(p,"-")·f)          uses invTᵀ base, no transpose
#   ( trp, "-")  D_p(invT_SE3(p,"-")ᵀ·f)         uses invTᵀ base, transposed
#
# NOTE  In cases (!trp,"+") and (!trp,"-"), the derivative of the upper block
# w.r.t. u contains a term involving T and tilde that is kept in comment and
# omitted from D_uu.  Only the φ-derivative part is retained (D_uu = D_φφ).
#
# EXTERNAL DEPENDENCIES
#   T_SE3_data, tilde, MAT3, MAT6, VEC3
#
# ==============================================================================

"""
    DinvT_SE3(a::T_SE3_data, f::VEC6; trp::Bool = false, sign_p::String = "+") → MAT6

Compute the 6×6 matrix of the Gâteaux derivative of `p ↦ invT_SE3(p)·f`
(or `p ↦ invT_SE3(p)ᵀ·f` when `trp = true`) with respect to the parameter
vector `p = [u; φ]`, using precomputed data `a` from [`T_SE3_input_data`](@ref).

The 6-vector `f` is split as `f_u = f[1:3]`, `f_φ = f[4:6]`.

Four combinations are available via `(trp, sign_p)`:

| `trp`  | `sign_p` | Computes |
|--------|----------|----------|
| `false` | `"+"` | `D_p(invT_SE3(p)·f)`           — forward, no transpose |
| `true`  | `"+"` | `D_p(invT_SE3(p)ᵀ·f)`          — forward, transposed   |
| `false` | `"-"` | `D_p(invT_SE3(p,"-")·f)`       — uses `invTᵀ` base     |
| `true`  | `"-"` | `D_p(invT_SE3(p,"-")ᵀ·f)`      — uses `invTᵀ` base, transposed |

    DinvT_SE3(p::VEC6, f::VEC6) → MAT6

Finite-difference approximation of `D_p(invT_SE3(p)·f)` (numerical verification).
"""
function DinvT_SE3(a::T_SE3_data, f::VEC6; trp::Bool = false, sign_p::String = "+")

    # Column-product helpers:
    #   pr3d(M, v)[i] = dM[i]·v   →  MAT3 with j-th column = M[j]·v
    #   pr3dt(M, v)   = same with each M[i] transposed
    pr3d(M, v)  = MAT3(M[1]*v,              M[2]*v,              M[3]*v)
    pr3dt(M, v) = MAT3(transpose(M[1])*v,   transpose(M[2])*v,   transpose(M[3])*v)

    invT  = a.invT
    dTu   = a.dTu
    Tu    = a.Tu
    dT    = a.dT
    dinvT = a.dinvT
    d²Tu  = a.d²Tu

    f_φ = VEC3(f[4], f[5], f[6])
    f_u = VEC3(f[1], f[2], f[3])

    # ── Case 1: D_p(invT_SE3(p)·f)  ──────────────────────────────────────────
    if !trp && sign_p == "+"

        g_φ      = invT * f_φ                           # T⁻¹(φ)·f_φ
        r        = -dTu * g_φ + tilde(Tu) * f_φ         # B·f_φ = invT·r
        dinvTf_φ = pr3d(dinvT, f_φ)                     # (∂T⁻¹/∂φᵢ)·f_φ columns  [3×3]
        d²Tug    = d²Tu[1]*g_φ[1] + d²Tu[2]*g_φ[2] + d²Tu[3]*g_φ[3]
        Dr_φ     = -tilde(f_φ) * dTu - d²Tug - dTu * dinvTf_φ

        D_φu = MAT3()
        D_φφ = dinvTf_φ
        # D_uu = D_φφ − invT·(dTg + tilde(f_φ)·T)  (tilde term omitted, see header note)
        D_uu = D_φφ
        D_uφ = pr3d(dinvT, f_u + r) + invT * Dr_φ

        return MAT6(D_uu, D_uφ, D_φu, D_φφ)

    # ── Case 2: D_p(invT_SE3(p)ᵀ·f)  ────────────────────────────────────────
    elseif trp && sign_p == "+"

        g_u   = transpose(invT) * f_u
        d²Tug = pr3dt(d²Tu, g_u)

        D_uu = MAT3()
        D_uφ = pr3dt(dinvT, f_u)
        D_φu = D_uφ
        D_φφ = -(tilde(Tu) + transpose(dTu * invT)) * pr3dt(dinvT, f_u) +
               pr3dt(dinvT, f_φ) + tilde(g_u) * dTu -
               pr3dt(dinvT, transpose(dTu) * g_u) - d²Tug * invT

        return MAT6(D_uu, D_uφ, D_φu, D_φφ)

    # ── Case 3: D_p(invT_SE3(p,"-")·f)  ─────────────────────────────────────
    elseif !trp && sign_p == "-"

        g_φ   = transpose(invT) * f_φ                   # (T⁻¹)ᵀ(φ)·f_φ
        d²Tug = d²Tu[1]*g_φ[1] + d²Tu[2]*g_φ[2] + d²Tu[3]*g_φ[3]

        D_φu = MAT3()
        D_φφ = pr3dt(dinvT, f_φ)
        # D_uu = D_φφ − (T⁻¹)ᵀ·dTg  (tilde term omitted, see header note)
        D_uu = D_φφ
        D_uφ = pr3dt(dinvT, f_u) -
               pr3dt(dinvT, dTu * g_φ) -
               transpose(invT) * (dTu * pr3dt(dinvT, f_φ) + d²Tug)

        return MAT6(D_uu, D_uφ, D_φu, D_φφ)

    # ── Case 4: D_p(invT_SE3(p,"-")ᵀ·f)  ────────────────────────────────────
    elseif trp && sign_p == "-"

        g_u = transpose(invT) * f_u

        D_uu = MAT3()
        D_uφ = pr3d(dinvT, f_u)
        D_φu = D_uφ
        D_φφ = pr3d(dinvT, f_φ) -
               pr3d(dinvT, transpose(dTu) * g_u) -
               invT * transpose(dTu) * pr3dt(dinvT, f_u) -
               invT * pr3dt(d²Tu, g_u)

        return MAT6(D_uu, D_uφ, D_φu, D_φφ)

    end
end

# Finite-difference validation: D_p(invT_SE3(p)·f) via centred differences
function DinvT_SE3(p::VEC6, f::VEC6)
    prec   = 1.0e-6
    pos    = VEC3(p[1], p[2], p[3])
    norm_x = norm2(pos)
    norm_x <= prec && (norm_x = 1.0)

    D = zeros(6, 6)
    for k = 1:6
        inc    = k <= 3 ? prec * norm_x : prec
        dp     = VEC6_unit(k, inc)
        col    = (0.5/inc) * (invT_SE3(p + dp) - invT_SE3(p - dp)) * f
        for i = 1:6
            D[i, k] = col[i]
        end
    end
    return MAT6(D)
end
