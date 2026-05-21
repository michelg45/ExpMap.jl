# ==============================================================================
# frame_interpol_2D.jl — SE(3) frame interpolation on a 2D surface element
# ==============================================================================
#
# TWO METHODS
#
# Method 1 — Reference surface interpolation
#   frame_interpol_2D(H, F, DF)
#   Interpolates the reference surface from nodal frames.  No precomputed
#   Jacobian is required: the surface Jacobian J_P is computed and returned.
#   Outputs: interpolated frame, physical strain vectors, tangent vectors,
#            unit surface normal, and the 2×2 surface Jacobian.
#
# Method 2 — Current (instantaneous) surface with sensitivity
#   frame_interpol_2D(H, F, DF, J)
#   Interpolates the current configuration using the reference Jacobian J
#   (returned by Method 1 for the same element) to map shape function
#   derivatives to physical coordinates.
#   Outputs: interpolated frame, physical strain vectors (displacement gradients
#            and curvatures), tangent vectors, unit normal, and the 2×Nnodes
#            array of nodal sensitivity matrices QN[α,k] = ∂f_P[α]/∂p_k.
#
# COMMON ALGORITHM (SE(3) weighted mean, Newton iterations)
#   Same as frame_interpol_1D, without the 2/ℓ scaling.
#   H_P = argmin Σ F[k] ‖log_SE3(H_P⁻¹·H[k])‖²
#   invA = Σ F[k]·invT_SE3(−p[k])
#
# STRAIN VECTORS
#   f_P[α] = invA⁻¹ · Σ_k DFs[k,α]·p[k]      α = 1, 2
#
#   Method 1: DFs = DF (reference; J_P is computed from f_P, then f_P is
#             transformed to physical coordinates via invJ = inv(J_P))
#   Method 2: DFs = DF·J⁻¹ (physical; DF is pre-scaled by the reference J)
#
# TANGENT VECTORS AND NORMAL
#   g_P[α] = R(H_P.phi) · VEC3(f_P[α][1:3])   (rotate translational part)
#   n_P    = (g_P[1] × g_P[2]) / (‖g_P[1]‖·‖g_P[2]‖)
#
# SENSITIVITY (Method 2 only)
#   QN[α,k] = ∂f_P[α]/∂p_k  (6×6 matrix)
#   See frame_interpol_1D for the chain-rule derivation.
#
# INPUT
#   H   :: Vector{NodeFrame}    nodal frames, length Nnodes
#   F   :: Vector{Float64}      shape function values at ξ
#   DF  :: Matrix{Float64}      shape function derivatives, size Nnodes×2
#   J   :: Matrix{Float64}      reference surface Jacobian, size 2×2 (Method 2)
#
# OUTPUT — Method 1
#   H_P :: NodeFrame
#   f_P :: Vector{VEC6}         physical strain vectors, length 2
#   g_P :: Vector{VEC3}         tangent vectors g[1], g[2], length 2
#   n_P :: VEC3                 unit surface normal
#   J_P :: Matrix{Float64}      2×2 surface Jacobian
#
# OUTPUT — Method 2
#   H_P :: NodeFrame
#   f_P :: Vector{VEC6}         physical strain vectors, length 2
#   g_P :: Vector{VEC3}         tangent vectors g[1], g[2], length 2
#   n_P :: VEC3                 unit surface normal
#   QN  :: Matrix{MAT6}         nodal sensitivities, size 2×Nnodes
#
# EXTERNAL DEPENDENCIES
#   NodeFrame, VEC3, VEC6, MAT6, log_SE3, exp_SE3, invT_SE3, R_SO3,
#   T_SE3_input_data, DinvT_SE3, crossp, norm2, zero, one
#
# ==============================================================================

"""
    frame_interpol_2D(H, F, DF) → (H_P, f_P, g_P, n_P, J_P)

Interpolate the **reference surface** from nodal frames `H[1:Nnodes]` at a
point with shape function values `F` and derivative matrix `DF` (size Nnodes×2).

Returns:
- `H_P :: NodeFrame`         — interpolated frame
- `f_P :: Vector{VEC6}`      — physical strain vectors in surface directions 1, 2
- `g_P :: Vector{VEC3}`      — tangent vectors `g[α] = R(φ)·u[α]`
- `n_P :: VEC3`              — unit surface normal `(g[1]×g[2])/‖…‖`
- `J_P :: Matrix{Float64}`   — 2×2 surface Jacobian `J_P[i,α] = f_P_ref[α][i]`

    frame_interpol_2D(H, F, DF, J) → (H_P, f_P, g_P, n_P, QN)

Interpolate the **current (instantaneous) surface** using the reference
Jacobian `J` (2×2, from the reference-configuration call above).  Shape
function derivatives are pre-scaled by `J⁻¹` before the strain assembly.

Returns additionally:
- `QN :: Matrix{MAT6}` — nodal sensitivities `QN[α,k] = ∂f_P[α]/∂p_k`,
  size 2×Nnodes.
"""
function frame_interpol_2D(H::Vector{NodeFrame}, F::Vector{Float64},
                            DF::Matrix{Float64})

    Nnodes = length(F)

    # ── Newton iteration for the SE(3) weighted mean ──────────────────────────
    H_P   = H[1]
    p     = Vector{VEC6}(undef, Nnodes)
    tol   = 1.0e-12
    test  = 2.0 * tol
    niter = 0
    nitmax = 20
    invA  = zero(MAT6)

    while test > tol && niter < nitmax

        r    = zero(VEC6)
        invA = zero(MAT6)

        for k = 1:Nnodes
            p[k] = log_SE3(-H_P * H[k])
            r    = r + F[k] * p[k]
        end

        test = norm2(r)
        niter == 0 && (tol = tol / test)

        for k = 1:Nnodes
            invA = invA + F[k] * invT_SE3(-p[k])
        end

        H_P   = H_P * exp_SE3(invA \ r)
        niter += 1
    end

    niter == nitmax && println("frame_interpol_2D: maximum iterations reached")

    # ── Strain vectors in reference coordinates ───────────────────────────────
    f_P_ref = Vector{VEC6}(undef, 2)
    f_P_ref[1] = zero(VEC6)
    f_P_ref[2] = zero(VEC6)
    for k = 1:Nnodes
        f_P_ref[1] = f_P_ref[1] + DF[k,1] * p[k]
        f_P_ref[2] = f_P_ref[2] + DF[k,2] * p[k]
    end
    f_P_ref[1] = invA \ f_P_ref[1]
    f_P_ref[2] = invA \ f_P_ref[2]

    # ── Reference surface Jacobian J_P and its inverse ────────────────────────
    # J_P[i,α] = i-th translational component of f_P_ref[α]  (i,α = 1,2)
    J_P = [f_P_ref[1][1]  f_P_ref[2][1];
           f_P_ref[1][2]  f_P_ref[2][2]]
    invJ = inv(J_P)

    # ── Transform strain vectors to physical surface coordinates ──────────────
    # f_P[α] = Σ_β  f_P_ref[β] · invJ[β,α]
    f_P = Vector{VEC6}(undef, 2)
    f_P[1] = f_P_ref[1]*invJ[1,1] + f_P_ref[2]*invJ[2,1]
    f_P[2] = f_P_ref[1]*invJ[1,2] + f_P_ref[2]*invJ[2,2]

    # ── Tangent vectors and unit normal ───────────────────────────────────────
    g_P = Vector{VEC3}(undef, 2)
    g_P[1] = R_SO3(H_P.phi, VEC3(f_P[1][1], f_P[1][2], f_P[1][3]))
    g_P[2] = R_SO3(H_P.phi, VEC3(f_P[2][1], f_P[2][2], f_P[2][3]))
    n_P    = crossp(g_P[1], g_P[2]) / (norm2(g_P[1]) * norm2(g_P[2]))

    return H_P, f_P, g_P, n_P, J_P
end


function frame_interpol_2D(H::Vector{NodeFrame}, F::Vector{Float64},
                            DF::Matrix{Float64}, J::Matrix{Float64})

    Nnodes = length(F)

    # Scale shape function derivatives from reference to physical coordinates
    # using the precomputed reference-surface Jacobian
    DFs = DF * inv(J)     # size Nnodes×2

    # ── Newton iteration for the SE(3) weighted mean ──────────────────────────
    H_P   = H[1]
    p     = Vector{VEC6}(undef, Nnodes)
    tol   = 1.0e-12
    test  = 2.0 * tol
    niter = 0
    nitmax = 20
    invA  = zero(MAT6)

    while test > tol && niter < nitmax

        r    = zero(VEC6)
        invA = zero(MAT6)

        for k = 1:Nnodes
            p[k] = log_SE3(-H_P * H[k])
            r    = r + F[k] * p[k]
        end

        test = norm2(r)
        niter == 0 && (tol = tol / test)

        for k = 1:Nnodes
            invA = invA + F[k] * invT_SE3(-p[k])
        end

        H_P   = H_P * exp_SE3(invA \ r)
        niter += 1
    end

    niter == nitmax && println("frame_interpol_2D: maximum iterations reached")

    # ── Physical strain vectors (DFs already includes J⁻¹ factor) ────────────
    f_P = Vector{VEC6}(undef, 2)
    f_P[1] = zero(VEC6)
    f_P[2] = zero(VEC6)
    for k = 1:Nnodes
        f_P[1] = f_P[1] + DFs[k,1] * p[k]
        f_P[2] = f_P[2] + DFs[k,2] * p[k]
    end
    f_P[1] = invA \ f_P[1]
    f_P[2] = invA \ f_P[2]

    # ── Tangent vectors and unit normal ───────────────────────────────────────
    g_P = Vector{VEC3}(undef, 2)
    g_P[1] = R_SO3(H_P.phi, VEC3(f_P[1][1], f_P[1][2], f_P[1][3]))
    g_P[2] = R_SO3(H_P.phi, VEC3(f_P[2][1], f_P[2][2], f_P[2][3]))
    n_P    = crossp(g_P[1], g_P[2]) / (norm2(g_P[1]) * norm2(g_P[2]))

    # ── Sensitivity matrices QN[α,k] = ∂f_P[α]/∂p_k ─────────────────────────
    # For each node k:
    #   invTm    = invT_SE3(p[k]; trp=true)     — transposed variant
    #   invTp    = invT_SE3(p[k]; trp=false)    — standard variant
    #   DTmpf[α] = D_p(invT_SE3⁻(p[k])·f_P[α]) — Gâteaux derivative, sign "-"
    #   Q[k]     = F[k]·invTm
    #   B[α]     = F[k]·DTmpf[α] + DFs[k,α]·I₆
    #   QN[α,k]  = B[α]·invTp
    #   C[α]    -= B[α]·invTm
    # Then: Q[k]    ← invA⁻¹·Q[k]
    #       QN[α,k] ← invA⁻¹·(QN[α,k] + C[α]·Q[k])

    I6  = one(MAT6)
    Q   = Vector{MAT6}(undef, Nnodes)
    QN  = Matrix{MAT6}(undef, 2, Nnodes)
    C   = [zero(MAT6), zero(MAT6)]

    for k = 1:Nnodes
        a      = T_SE3_input_data(p[k])
        invTm  = invT_SE3(a; trp = true)
        invTp  = invT_SE3(a; trp = false)
        DTmpf1 = DinvT_SE3(a, f_P[1]; sign_p = "-")
        DTmpf2 = DinvT_SE3(a, f_P[2]; sign_p = "-")
        Q[k]   = F[k] * invTm
        B1     = F[k] * DTmpf1 + DFs[k,1] * I6
        B2     = F[k] * DTmpf2 + DFs[k,2] * I6
        QN[1,k] = B1 * invTp
        QN[2,k] = B2 * invTp
        C[1]   = C[1] - B1 * invTm
        C[2]   = C[2] - B2 * invTm
    end

    for k = 1:Nnodes
        Q[k]    = invA \ Q[k]
        QN[1,k] = invA \ (QN[1,k] + C[1] * Q[k])
        QN[2,k] = invA \ (QN[2,k] + C[2] * Q[k])
    end

    return H_P, f_P, g_P, n_P, QN
end
