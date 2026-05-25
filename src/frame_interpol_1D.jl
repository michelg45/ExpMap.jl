# ==============================================================================
# frame_interpol_1D.jl — SE(3) frame interpolation on a 1D element
# ==============================================================================
#
# Computes the interpolated NodeFrame H_P at a given point ξ of a 1D element
# using the SE(3)-weighted mean of the nodal frames, then derives the
# generalized strain vector f_P and the nodal sensitivity matrices QN[k].
#
# ALGORITHM
#   The interpolated frame minimises the weighted sum of squared relative
#   motions in the Lie algebra sense.  It is found by Newton iteration on SE(3):
#
#     1. Initialise H_P = H[1]
#     2. For each node k: p[k] = log_SE3(H_P⁻¹ · H[k])
#     3. Residual:  r = Σ F[k] · p[k]
#     4. Tangent:   A = Σ F[k] · invT_SE3(−p[k])
#     5. Update:    H_P ← H_P · exp_SE3(A⁻¹·r)
#     6. Repeat until ‖r‖ < ε_rel · ‖r₀‖
#
# STRAIN VECTOR
#   After convergence, the generalised strain vector is:
#
#     f_P = A⁻¹ · Σ DF[k] · p[k]
#
#   where DF[k] = (2/ℓ) · (dN_k/dξ) are the shape function derivatives
#   mapped from the reference segment [−1,+1] to the physical length ℓ.
#   Components 1–3 are displacement-gradient measures; 4–6 are curvatures.
#
# SENSITIVITY MATRICES
#   QN[k] = ∂f_P/∂p_k   (6×6, one per node)
#
#   computed via the chain rule using the Gâteaux derivatives DinvT_SE3.
#
# INPUT
#   H  :: Vector{NodeFrame}   nodal frames, length Nnodes
#   F  :: Vector{Float64}     Lagrange shape function values at ξ
#   DF :: Vector{Float64}     shape function derivatives dN_k/dξ at ξ
#   ℓ  :: Float64             element length (physical, for the 2/ℓ scaling)
#
# OUTPUT
#   H_P :: NodeFrame           interpolated frame
#   f_P :: VEC6                generalised strain vector
#   QN  :: Vector{MAT6}        sensitivity matrices ∂f_P/∂p_k, length Nnodes
#
# EXTERNAL DEPENDENCIES
#   NodeFrame, VEC6, MAT6, log_SE3, exp_SE3, invT_SE3, T_SE3_input_data,
#   DinvT_SE3, norm2, zero, one
#
# ==============================================================================

"""
    frame_interpol_1D(H, F, DF, ℓ) → (H_P, f_P, QN)

Interpolate a frame curve from nodal frames `H[1:Nnodes]` at a point with
shape function values `F` and derivatives `DF` on a 1D element of length `ℓ`.

Returns:
- `H_P :: NodeFrame`    — interpolated frame (SE(3) weighted mean)
- `f_P :: VEC6`         — generalised strains: `[displacement gradients; curvatures]`
- `QN  :: Vector{MAT6}` — nodal sensitivity matrices `∂f_P/∂p_k`

The shape function derivatives `DF` are given in reference coordinates ξ ∈ [−1,+1]
and scaled internally by `2/ℓ` to obtain physical derivatives `d/ds`.
"""
function frame_interpol_1D(H::Vector{NodeFrame}, F::Vector{Float64},
                            DF::Vector{Float64}, ℓ::Float64)

    Nnodes = length(F)

    # Shape function derivatives scaled to physical arc length
    # ξ ∈ [−1,+1]  →  s ∈ [0,ℓ]:  d/ds = (2/ℓ)·d/dξ
    DFs = DF .* (2.0/ℓ)

    # ── Newton iteration for the SE(3) weighted mean ──────────────────────────
    H_P   = H[1]          # initial guess
    p     = Vector{VEC6}(undef, Nnodes)
    tol   = 1.0e-12
    test  = 2.0 * tol     # dummy value to enter the loop
    niter = 0
    nitmax = 20
    invA  = zero(MAT6)

    while test > tol && niter < nitmax

        r    = zero(VEC6)
        invA = zero(MAT6)

        for k = 1:Nnodes
            p[k] = log_SE3(-H_P * H[k])          # relative motion in se(3)
            r    = r + F[k] * p[k]               # weighted residual
        end

        test = norm2(r)
        niter == 0 && (tol = tol / test)          # switch to relative tolerance

        for k = 1:Nnodes
            invA = invA + F[k] * invT_SE3(-p[k]) # tangent accumulation
        end

        H_P   = H_P * exp_SE3(invA \ r)          # SE(3) update
        niter += 1
    end

    niter == nitmax && println("frame_interpol_1D: maximum iterations reached")

    # ── Generalised strain vector ─────────────────────────────────────────────
    f_P = zero(VEC6)
    for k = 1:Nnodes
        f_P = f_P + DFs[k] * p[k]
    end
    f_P = invA \ f_P

    # ── Sensitivity matrices Q[k] = ∂p/∂p_k  and  QN[k] = ∂f_P/∂p_k ─────────
    # For each node k (at converged p[k]):
    #   a     = T_SE3_input_data(p[k])         — precomputed SO(3) data
    #   invTm = invT_SE3(a, "-")               — invT_SE3(−p[k])
    #   invTp = invT_SE3(a, "+")               — invT_SE3(p[k])
    #   DTmpf = D_p(invT_SE3(−p[k])·f_P)      — Gâteaux derivative at −p
    #   Q[k]  = F[k]·invTp
    #   B     = F[k]·DTmpf + DFs[k]·I₆
    #   QN[k] = B·invTp                        (before left-multiplication by A⁻¹)
    #   C    -= B·invTm                        (accumulated correction)
    # Then: Q[k]  ← A⁻¹·Q[k]
    #       QN[k] ← A⁻¹·(QN[k] + C·Q[k])

    I6 = one(MAT6)

    Q  = Vector{MAT6}(undef, Nnodes)
    QN = Vector{MAT6}(undef, Nnodes)
    C  = zero(MAT6)

    for k = 1:Nnodes
        a     = T_SE3_input_data(p[k])
        invTm = invT_SE3(a, "-")                        # invT_SE3(−p[k])
        invTp = invT_SE3(a, "+")                        # invT_SE3(p[k])
        DTmpf = DinvT_SE3(a, f_P; sign_p = "-")        # D_p(invT_SE3(−p[k])·f_P)
        Q[k]  = F[k] * invTp
        B     = F[k] * DTmpf + DFs[k] * I6
        QN[k] = B * invTp
        C     = C - B * invTm
    end

    for k = 1:Nnodes
        Q[k]  = invA \ Q[k]
        QN[k] = invA \ (QN[k] + C * Q[k])
    end

    return H_P, f_P, QN
end
