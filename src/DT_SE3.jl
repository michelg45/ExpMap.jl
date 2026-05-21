# ==============================================================================
# DT_SE3.jl — Gâteaux Derivative of the SE(3) Tangent Operator
# ==============================================================================
#
# Computes the 6×6 matrix of the Gâteaux derivative of  p ↦ T_SE3(p)·f
# (or  p ↦ T_SE3(p)ᵀ·f)  with respect to the full parameter vector p = [u; φ].
#
# INPUT
#   a       :: T_SE3_data   precomputed data from T_SE3_input_data
#   f       :: VEC6         fixed 6-vector; split as f_u = f[1:3], f_φ = f[4:6]
#   trp     :: Bool = false transpose flag (true → differentiate T_SE3ᵀ·f)
#   sign_p  :: String = "+" operator sign ("+" → forward T_SE3; "-" → uses Tᵀ)
#
# OUTPUT   MAT6   — the 6×6 Gâteaux derivative matrix written in block form:
#
#   [ D_uu  D_uφ ]     rows 1–3: derivative of upper-block output
#   [ D_φu  D_φφ ]     rows 4–6: derivative of lower-block output
#
# FOUR CASES   (trp) × (sign_p)
#
#   (!trp, "+")  D_p(T_SE3(p)·f)              forward, no transpose
#   ( trp, "+")  D_p(T_SE3(p)ᵀ·f)            forward, transposed
#   (!trp, "-")  D_p(T_SE3(p,"-")·f)          uses Tᵀ base, no transpose
#   ( trp, "-")  D_p(T_SE3(p,"-")ᵀ·f)         uses Tᵀ base, transposed
#
# NOTE  In cases (!trp,"+") and (!trp,"-"), the derivative of the upper block
# w.r.t. u contains a term  tilde(T·f_φ)·T  that is kept in comment and
# omitted from D_uu.  Only the φ-derivative part is retained.
#
# EXTERNAL DEPENDENCIES
#   T_SE3_data, tilde, MAT3, MAT6, VEC3
#
# ==============================================================================

"""
    DT_SE3(a::T_SE3_data, f::VEC6; trp::Bool = false, sign_p::String = "+") → MAT6

Compute the 6×6 matrix of the Gâteaux derivative of `p ↦ T_SE3(p)·f`
(or `p ↦ T_SE3(p)ᵀ·f` when `trp = true`) with respect to the parameter
vector `p = [u; φ]`, using precomputed data `a` from [`T_SE3_input_data`](@ref).

The 6-vector `f` is split as `f_u = f[1:3]`, `f_φ = f[4:6]`.

Four combinations are available via `(trp, sign_p)`:

| `trp`  | `sign_p` | Computes |
|--------|----------|----------|
| `false` | `"+"` | `D_p(T_SE3(p)·f)`           — forward, no transpose |
| `true`  | `"+"` | `D_p(T_SE3(p)ᵀ·f)`          — forward, transposed   |
| `false` | `"-"` | `D_p(T_SE3(p,"-")·f)`       — uses `Tᵀ` base       |
| `true`  | `"-"` | `D_p(T_SE3(p,"-")ᵀ·f)`      — uses `Tᵀ` base, transposed |
"""
function DT_SE3(a::T_SE3_data, f::VEC6; trp::Bool = false, sign_p::String = "+")

    # Column-product helpers:
    #   pr3d(M, v)[i] = dM[i]·v   →  MAT3 with j-th column = M[j]·v
    #   pr3dt(M, v)   = same with each M[i] transposed
    pr3d(M, v)  = MAT3(M[1]*v,              M[2]*v,              M[3]*v)
    pr3dt(M, v) = MAT3(transpose(M[1])*v,   transpose(M[2])*v,   transpose(M[3])*v)

    T    = a.T
    Tu   = a.Tu
    dTu  = a.dTu
    dT   = a.dT
    u    = VEC3(a.p[1], a.p[2], a.p[3])

    f_φ = VEC3(f[4], f[5], f[6])
    f_u = VEC3(f[1], f[2], f[3])

    # ── Case 1: D_p(T_SE3(p)·f)  ─────────────────────────────────────────────
    if !trp && sign_p == "+"

        g_φ   = T * f_φ                             # T(φ)·f_φ
        dTf_φ = pr3d(dT, f_φ)                       # (∂T/∂φᵢ)·f_φ columns  [3×3]

        D_φu = MAT3()
        D_φφ = dTf_φ
        # D_uu = dTf_φ + tilde(g_φ)*T  (tilde term omitted, see header note)
        D_uu = D_φφ
        D_uφ = pr3d(dT, f_u) + pr3d(a.d²Tu, f_φ) - tilde(Tu)*dTf_φ + tilde(g_φ)*dTu

        return MAT6(D_uu, D_uφ, D_φu, D_φφ)

    # ── Case 2: D_p(T_SE3(p)ᵀ·f)  ───────────────────────────────────────────
    elseif trp && sign_p == "+"

        D_uu = MAT3()
        D_uφ = pr3dt(dT, f_u)
        D_φu = D_uφ
        D_φφ = pr3dt(a.d²Tu, f_u) + pr3dt(dT, f_φ) -
               transpose(T) * tilde(f_u) * dTu - pr3dt(dT, tilde(f_u)*Tu)

        return MAT6(D_uu, D_uφ, D_φu, D_φφ)

    # ── Case 3: D_p(T_SE3(p,"-")·f)  ────────────────────────────────────────
    elseif !trp && sign_p == "-"

        g_φ   = transpose(T) * f_φ                  # Tᵀ(φ)·f_φ
        dTf_φ = pr3dt(dT, f_φ)                      # (∂T/∂φᵢ)ᵀ·f_φ columns
        dTtu  = pr3dt(dT, u)                         # (∂T/∂φᵢ)ᵀ·u  columns

        D_φu = MAT3()
        D_φφ = dTf_φ
        # D_uu = dTf_φ - tilde(g_φ)*Tᵀ  (tilde term omitted, see header note)
        D_uu = D_φφ
        D_uφ = pr3dt(dT, f_u) + pr3d(a.d²Ttu, f_φ) +
               tilde(transpose(T)*u) * dTf_φ - tilde(g_φ) * dTtu

        return MAT6(D_uu, D_uφ, D_φu, D_φφ)

    # ── Case 4: D_p(T_SE3(p,"-")ᵀ·f)  ──────────────────────────────────────
    elseif trp && sign_p == "-"

        D_uu = MAT3()
        D_uφ = pr3d(dT, f_u)
        D_φu = D_uφ
        D_φφ = pr3d(dT, f_φ) +
               pr3d(dT, tilde(f_u) * transpose(T) * u) +
               T * tilde(f_u) * pr3dt(dT, u) +
               pr3dt(a.d²Ttu, f_u)

        return MAT6(D_uu, D_uφ, D_φu, D_φφ)

    end
end
