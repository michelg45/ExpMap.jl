# ==============================================================================
# T_SE3_input_data.jl — Complete data for the SE(3) tangent operator
# ==============================================================================
#
# Computes all quantities stored in T_SE3_data: the SO(3) tangent operator T(φ)
# and its inverse T⁻¹(φ), their first and second derivatives w.r.t. φ, and the
# derived vectors/matrices involving the translation part u.
#
# INPUT
#   p::VEC6 — 6-vector [u; φ], with translation u = p[1:3] and
#             Cartesian rotation vector φ = p[4:6].
#
# OUTPUT   T_SE3_data struct with fields:
#   p       :: VEC6             input vector
#   invT    :: MAT3             T⁻¹(φ) = I + ½φ̃ + γ(x)·W
#   T       :: MAT3             T(φ)   = I + α(x)·tilde(φ̂) + β(x)·W
#   Tu      :: VEC3             T(φ)·u
#   dTu     :: MAT3             Gâteaux derivative matrix of  φ ↦ T(φ)·u
#   dinvT   :: Vector{MAT3}[3]  dinvT[i] = ∂T⁻¹/∂φᵢ
#   dT      :: Vector{MAT3}[3]  dT[i]    = ∂T/∂φᵢ
#   d²Tu    :: Vector{MAT3}[3]  d²Tu[i][:,j]  = (∂²T/∂φᵢ∂φⱼ)·u
#   d²Ttu   :: Vector{MAT3}[3]  d²Ttu[i][:,j] = (∂²Tᵀ/∂φᵢ∂φⱼ)·u
#
# NOTATION (consistent with T_SO3, DT_SO3, invT_SO3)
#   x  = ‖φ‖,   φ̂ = φ/x,   W = φ̂⊗φ̂ − I
#   α, β, γ    scalar coefficients from T_functions, invT_functions
#   dφ̂ = (I − φ̂⊗φ̂)/x = −W/x          Jacobian of (φ → φ̂)       [3×3]
#   dφ̂[i,j]                            Hessian of ‖φ‖  (same matrix)
#   dW[i]  = ∂W/∂φᵢ  = φ̂⊗dφ̂[:,i] + dφ̂[:,i]⊗φ̂              [3×3]
#   d²φ̂[k][:,j] = ∂²φ̂/∂φₖ∂φⱼ
#              = (−φ̂ₖI − eₖ⊗φ̂ − φ̂⊗eₖ + 3φ̂ₖ·φ̂⊗φ̂) / x²
#   d²W[i,j] = ∂²W/∂φᵢ∂φⱼ
#
# SMALL-ANGLE LIMIT  ‖φ‖² < 1e-8
#   T ≈ I−½φ̃,  T⁻¹ ≈ I+½φ̃,  dT[i] = −½tilde(eᵢ),  d² terms → 0
#
# EXTERNAL DEPENDENCIES
#   T_functions(x; der=2), invT_functions(x; der=1)
#   tilde, outerp, crossp, dotp, one(MAT3), MAT3, VEC3
#
# ==============================================================================

"""
    T_SE3_input_data(p::VEC6) → T_SE3_data

Compute all quantities in [`T_SE3_data`](@ref) needed to assemble the SE(3)
tangent operator and its directional derivative `D_p(T_SE3(p)·a)`.

Given `p = [u; φ]` with translation `u = p[1:3]` and rotation vector
`φ = p[4:6]`, returns a `T_SE3_data` containing:

- `T`, `invT`           — SO(3) tangent operator and its inverse
- `Tu`                  — `T(φ)·u`
- `dTu`                 — Gâteaux derivative matrix of `φ ↦ T(φ)·u`
- `dT[i]`, `dinvT[i]`   — partial derivatives `∂T/∂φᵢ`, `∂T⁻¹/∂φᵢ`
- `d²Tu[i][:,j]`        — `(∂²T/∂φᵢ∂φⱼ)·u`
- `d²Ttu[i][:,j]`       — `(∂²Tᵀ/∂φᵢ∂φⱼ)·u`

All shared intermediate quantities (`φ̂`, `α`, `β`, `γ` and their derivatives,
`W`, `dφ̂`, `dW`, `d²φ̂`, `d²W`) are computed once and reused throughout.
"""
function T_SE3_input_data(p::VEC6)
    φ  = VEC3(p[4], p[5], p[6])   # rotation vector
    u  = VEC3(p[1], p[2], p[3])   # translation vector

    I  = one(MAT3)
    x² = dotp(φ, φ)

    dT    = Vector{MAT3}(undef, 3)
    dinvT = Vector{MAT3}(undef, 3)
    d²Tu  = Vector{MAT3}(undef, 3)
    d²Ttu = Vector{MAT3}(undef, 3)

    # Unit vectors e₁, e₂, e₃ (columns of I) — for tilde(eᵢ) and outer products
    e = (I[:,1], I[:,2], I[:,3])

    # ── Small-angle limit ─────────────────────────────────────────────────────
    if x² < 1.0e-8
        T    = I - 0.5 * tilde(φ)        # T(φ)   ≈ I − ½φ̃
        invT = I + 0.5 * tilde(φ)        # T⁻¹(φ) ≈ I + ½φ̃
        Tu   = u - 0.5 * crossp(φ, u)    # T(φ)·u ≈ u − ½(φ×u)
        dTu  = 0.5 * tilde(u)            # Gâteaux derivative ≈ ½ã
        for i = 1:3
            dT[i]    = -0.5 * tilde(e[i])
            dinvT[i] =  0.5 * tilde(e[i])
            d²Tu[i]  = MAT3()
            d²Ttu[i] = MAT3()
        end
        return T_SE3_data(p, invT, T, Tu, dTu, dinvT, dT, d²Tu, d²Ttu)
    end

    # ── Shared intermediates ──────────────────────────────────────────────────
    x   = sqrt(x²)
    φ̂   = VEC3(φ[1]/x, φ[2]/x, φ[3]/x)        # unit rotation axis
    # gradient of ‖φ‖: d‖φ‖/dφ = φ̂  (φ̂[i] used as φ̂ᵢ below)

    α, dα, ddα, β, dβ, ddβ = T_functions(x; der=2)
    γ, dγ                   = invT_functions(x; der=1)

    Opp = outerp(φ̂, φ̂)                         # φ̂⊗φ̂
    W   = Opp - I                               # φ̂⊗φ̂ − I
    Wu  = W * u                                 # (φ̂⊗φ̂ − I)·u
    φ̂xu = crossp(φ̂, u)                          # φ̂×u

    T    = I + α * tilde(φ̂) + β * W
    invT = I + 0.5 * tilde(φ) + γ * W
    Tu   = u + α * φ̂xu + β * Wu

    # Jacobian of (φ → φ̂):  dφ̂ = (I − φ̂⊗φ̂)/x = −W/x
    # Its (i,j) entry is also the Hessian of ‖φ‖: ∂²‖φ‖/∂φᵢ∂φⱼ = dφ̂[i,j]
    dφ̂ = (I - Opp) / x

    # ── First derivatives of W, T, T⁻¹ ──────────────────────────────────────
    dW = Vector{MAT3}(undef, 3)
    for i = 1:3
        dW[i] = outerp(φ̂, dφ̂[:,i]) + outerp(dφ̂[:,i], φ̂)    # ∂W/∂φᵢ
    end

    for i = 1:3
        φ̂ᵢ      = φ̂[i]                                       # = d‖φ‖/dφᵢ
        dT[i]    = α * tilde(dφ̂[:,i]) + β * dW[i] +
                   φ̂ᵢ * (dα * tilde(φ̂) + dβ * W)
        dinvT[i] = 0.5 * tilde(e[i]) + γ * dW[i] + φ̂ᵢ * dγ * W
    end

    dTu = MAT3(dT[1]*u, dT[2]*u, dT[3]*u)

    # ── Second derivatives ────────────────────────────────────────────────────
    # Hessian of φ̂:  d²φ̂[k][:,j] = ∂²φ̂/∂φₖ∂φⱼ
    #   = (−φ̂ₖI − eₖ⊗φ̂ − φ̂⊗eₖ + 3φ̂ₖ·(φ̂⊗φ̂)) / x²
    d²φ̂ = Vector{MAT3}(undef, 3)
    for k = 1:3
        φ̂ₖ = φ̂[k]
        d²φ̂[k] = (-(φ̂ₖ * I + outerp(e[k], φ̂) + outerp(φ̂, e[k])) + 3φ̂ₖ * Opp) / x²
    end

    # Second derivatives of W and T (lower triangle, then symmetrised)
    d²W = Array{MAT3}(undef, 3, 3)
    d²T = Array{MAT3}(undef, 3, 3)

    for i = 1:3, j = 1:i
        # ∂²W/∂φᵢ∂φⱼ = (∂φ̂/∂φⱼ)⊗(∂φ̂/∂φᵢ) + (∂φ̂/∂φᵢ)⊗(∂φ̂/∂φⱼ)
        #              + φ̂⊗(∂²φ̂/∂φᵢ∂φⱼ) + (∂²φ̂/∂φᵢ∂φⱼ)⊗φ̂
        d²φ̂ⱼᵢ  = d²φ̂[j][:,i]                 # ∂²φ̂/∂φᵢ∂φⱼ as a VEC3
        d²W[i,j] = outerp(dφ̂[:,j], dφ̂[:,i]) + outerp(dφ̂[:,i], dφ̂[:,j]) +
                   outerp(φ̂, d²φ̂ⱼᵢ) + outerp(d²φ̂ⱼᵢ, φ̂)

        φ̂ᵢ  = φ̂[i];  φ̂ⱼ = φ̂[j]
        Hxij = dφ̂[i,j]               # ∂²‖φ‖/∂φᵢ∂φⱼ = dφ̂[i,j]

        d²T[i,j] = α   * tilde(d²φ̂ⱼᵢ) +
                   dα  * (φ̂ⱼ * tilde(dφ̂[:,i]) + Hxij * tilde(φ̂) + φ̂ᵢ * tilde(dφ̂[:,j])) +
                   ddα * φ̂ᵢ * φ̂ⱼ * tilde(φ̂) +
                    β  * d²W[i,j] +
                   dβ  * (φ̂ⱼ * dW[i] + Hxij * W + φ̂ᵢ * dW[j]) +
                   ddβ * φ̂ᵢ * φ̂ⱼ * W

        if j != i
            d²W[j,i] = d²W[i,j]
            d²T[j,i] = d²T[i,j]
        end
    end

    # ── Build d²Tu and d²Ttu from d²T ────────────────────────────────────────
    for i = 1:3
        d²Tu[i]  = MAT3(d²T[i,1]*u, d²T[i,2]*u, d²T[i,3]*u)
        d²Ttu[i] = MAT3(transpose(d²T[i,1])*u,
                        transpose(d²T[i,2])*u,
                        transpose(d²T[i,3])*u)
    end

    return T_SE3_data(p, invT, T, Tu, dTu, dinvT, dT, d²Tu, d²Ttu)
end
