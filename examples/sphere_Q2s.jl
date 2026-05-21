# ==============================================================================
# sphere_Q2s.jl — SE(3) frame interpolation on a spherical Q2s element
# ==============================================================================
#
# Geometry:
#   Sphere of radius R = 10, centred at C = [2, 0, 0].
#   The element occupies the patch cut by  z ∈ [6, 8]  and  φ ∈ [0, π/2].
#
#   At z = 6:  (x−2)²+y² = 64  ⟹  ρ₁ = 8,  cos θ₁ = 3/5,  sin θ₁ = 4/5
#   At z = 8:  (x−2)²+y² = 36  ⟹  ρ₂ = 6,  cos θ₂ = 4/5,  sin θ₂ = 3/5
#   θ₁ + θ₂ = π/2  (3-4-5 Pythagorean triplet), so θ_m = π/4.
#
# Node numbering (Q2s, 8 nodes):
#
#         φ=0          φ=π/4        φ=π/2
#   z=8   H[4] ─────── H[7] ─────── H[3]
#          │                          │
#   z=mid  H[8]                      H[6]
#          │                          │
#   z=6   H[1] ─────── H[5] ─────── H[2]
#
#   Corners (θ, φ):
#     H[1]: (θ₁, 0)      H[2]: (θ₁, π/2)
#     H[3]: (θ₂, π/2)    H[4]: (θ₂, 0)
#   Mid-side nodes:
#     H[5]: (θ₁, π/4)  — mid of H[1]–H[2]
#     H[6]: (π/4, π/2) — mid of H[2]–H[3]   (θ_m = π/4)
#     H[7]: (θ₂, π/4)  — mid of H[3]–H[4]
#     H[8]: (π/4,  0)  — mid of H[4]–H[1]   (θ_m = π/4)
#
# Local frame at each node:
#   e_x = e_φ  =  (−sin φ,   cos φ,  0  )   horizontal tangent
#   e_y = −e_θ =  (−cos θ cos φ, −cos θ sin φ, sin θ)   meridional tangent
#   e_z = n̂   =  ( sin θ cos φ,  sin θ sin φ, cos θ)   outward normal
#
# The element reference domain for shape_functions_2D(8, ξ) is [−1,+1]² (Q2s).
# ==============================================================================

if !@isdefined(ExpMap)
    include(joinpath(@__DIR__, "..", "src", "ExpMap.jl"))
    using .ExpMap
end
using Printf
using LinearAlgebra: det

# ── Sphere parameters ─────────────────────────────────────────────────────────
const R_sphere  = 10.0
const C_sphere  = VEC3(2.0, 0.0, 0.0)

# ── Helper: NodeFrame for a point on the sphere at (θ, φ) ─────────────────────
function sphere_frame(θ::Float64, φ::Float64)
    sθ, cθ = sin(θ), cos(θ)
    sφ, cφ = sin(φ), cos(φ)

    # Position on sphere
    x = C_sphere + R_sphere * VEC3(sθ*cφ, sθ*sφ, cθ)

    # Local orthonormal frame
    e_x = VEC3(-sφ,       cφ,       0.0 )   # e_φ  — horizontal tangent
    e_y = VEC3(-cθ*cφ,  -cθ*sφ,    sθ  )   # −e_θ — meridional tangent
    e_z = VEC3( sθ*cφ,   sθ*sφ,    cθ  )   # n̂   — outward normal

    ψ = invR_SO3(MAT3(e_x, e_y, e_z))
    return NodeFrame(x, ψ)
end

# ── Polar angles for the two latitude circles ──────────────────────────────────
# z = 6: cos θ₁ = 6/10 = 3/5,  sin θ₁ = 8/10 = 4/5
# z = 8: cos θ₂ = 8/10 = 4/5,  sin θ₂ = 6/10 = 3/5
θ₁ = acos(3.0/5.0)   # ≈ 53.13°  (lower circle, z=6)
θ₂ = acos(4.0/5.0)   # ≈ 36.87°  (upper circle, z=8)
θ_m = π/4            # mid-latitude angle  (exact: θ₁+θ₂ = π/2)

# ── 8 nodal frames ────────────────────────────────────────────────────────────
H = Vector{NodeFrame}(undef, 8)

# Corners
H[1] = sphere_frame(θ₁,   0.0  )   # (θ₁, 0)
H[2] = sphere_frame(θ₁,   π/2  )   # (θ₁, π/2)
H[3] = sphere_frame(θ₂,   π/2  )   # (θ₂, π/2)
H[4] = sphere_frame(θ₂,   0.0  )   # (θ₂, 0)

# Mid-side nodes
H[5] = sphere_frame(θ₁,   π/4  )   # mid H[1]–H[2]
H[6] = sphere_frame(θ_m,  π/2  )   # mid H[2]–H[3]
H[7] = sphere_frame(θ₂,   π/4  )   # mid H[3]–H[4]
H[8] = sphere_frame(θ_m,  0.0  )   # mid H[4]–H[1]

# ── Helper: print a NodeFrame (position + rotation) ───────────────────────────
function print_frame(label::String, h::NodeFrame)
    x  = h.x
    ψ  = h.phi
    R  = R_SO3(ψ)           # rotation matrix (columns = local e_x, e_y, e_z)
    @printf("  %s\n", label)
    @printf("    x   = [%8.5f, %8.5f, %8.5f]\n", x[1],    x[2],    x[3])
    @printf("    ψ   = [%8.5f, %8.5f, %8.5f]   (RV3, ‖ψ‖ = %.5f)\n",
            ψ[1], ψ[2], ψ[3], norm2(VEC3(ψ[1],ψ[2],ψ[3])))
    @printf("    e_x = [%8.5f, %8.5f, %8.5f]   (horizontal tangent e_φ)\n",
            R[1,1], R[2,1], R[3,1])
    @printf("    e_y = [%8.5f, %8.5f, %8.5f]   (meridional tangent −e_θ)\n",
            R[1,2], R[2,2], R[3,2])
    @printf("    e_z = [%8.5f, %8.5f, %8.5f]   (outward normal n̂)\n",
            R[1,3], R[2,3], R[3,3])
end

# ── Print all nodal frames ─────────────────────────────────────────────────────
println("=== Q2s element on sphere (R=10, C=[2,0,0]) ===\n")
println("Nodal frames — position x, rotation vector ψ, and frame axes:\n")
print_frame("H[1]  (θ₁,  0  ) — corner",     H[1])
print_frame("H[2]  (θ₁,  π/2) — corner",     H[2])
print_frame("H[3]  (θ₂,  π/2) — corner",     H[3])
print_frame("H[4]  (θ₂,  0  ) — corner",     H[4])
print_frame("H[5]  (θ₁,  π/4) — mid H1–H2",  H[5])
print_frame("H[6]  (π/4, π/2) — mid H2–H3",  H[6])
print_frame("H[7]  (θ₂,  π/4) — mid H3–H4",  H[7])
print_frame("H[8]  (π/4, 0  ) — mid H4–H1",  H[8])

# ── Verification: all nodes should lie on the sphere ──────────────────────────
println("\nDistance to sphere surface (should all be ≈ 0):")
for (i, h) in enumerate(H)
    d = norm2(h.x - C_sphere) - R_sphere
    @printf("  H[%d]: ‖x−C‖ − R = %+.2e\n", i, d)
end

# ── Reference surface interpolation ───────────────────────────────────────────
println("\n=== Reference surface interpolation at centre point ξ = [0, 0] ===\n")

ξ_c = [0.0, 0.0]
F_c, DF_c = shape_functions_2D(8, ξ_c)
H_P_ref, f_P_ref, g_P_ref, n_P_ref, J_P = frame_interpol_2D(H, F_c, DF_c)

println("Interpolated frame at ξ = (0, 0):")
println("  position:  $(H_P_ref.x)")
println("  rot. vec.: $(H_P_ref.phi)")

println("\nTangent vectors g_P_ref:")
println("  g_P[1] = $(g_P_ref[1])")
println("  g_P[2] = $(g_P_ref[2])")
println("Unit normal n_P = $(n_P_ref)")

println("\nSurface Jacobian J_P:")
println("  $(J_P[1,1])  $(J_P[1,2])")
println("  $(J_P[2,1])  $(J_P[2,2])")

# ── 2×2 Gauss quadrature on the element ───────────────────────────────────────
println("\n=== 2×2 Gauss quadrature on current configuration ===\n")

ξ_gauss, w_gauss = gauss_points(2, 4)   # 4 = 2×2 points on [−1,+1]²
Ngauss = length(w_gauss)

for ig = 1:Ngauss
    F_g, DF_g = shape_functions_2D(8, ξ_gauss[ig])

    # Reference interpolation provides J for this element (same J for all
    # Gauss points since the reference configuration is fixed; recomputed
    # here for completeness — in production, use the J from one call)
    _, _, _, _, J_ig = frame_interpol_2D(H, F_g, DF_g)

    # Current configuration with sensitivities (H == reference here)
    H_P, f_P, g_P, n_P, QN = frame_interpol_2D(H, F_g, DF_g, J_ig)

    dA = det(J_ig) * w_gauss[ig]

    @printf("  Gauss point %d: ξ = [%.4f, %.4f],  w = %.4f,  dA = %.4f\n",
            ig, ξ_gauss[ig][1], ξ_gauss[ig][2], w_gauss[ig], dA)
    @printf("    position:  [%.4f, %.4f, %.4f]\n",
            H_P.x[1], H_P.x[2], H_P.x[3])
    @printf("    normal:    [%.4f, %.4f, %.4f]\n",
            n_P[1], n_P[2], n_P[3])
    @printf("    ‖n_P‖ − 1 = %+.2e\n", norm2(n_P) - 1.0)
end

# ── Total area estimate ────────────────────────────────────────────────────────
println("\n=== Total element area (numerical vs analytical) ===\n")

area_num = 0.0
for ig = 1:Ngauss
    global area_num
    F_g, DF_g = shape_functions_2D(8, ξ_gauss[ig])
    _, _, g_P, _, J_ig = frame_interpol_2D(H, F_g, DF_g)
    # Area element = ‖g_P[1] × g_P[2]‖ · |det J| · w
    area_num += norm2(crossp(g_P[1], g_P[2])) * abs(det(J_ig)) * w_gauss[ig]
end

# Analytical area of spherical patch: φ ∈ [0,π/2], θ ∈ [θ₂, θ₁]
# A = R² ∫₀^{π/2} dφ ∫_{θ₂}^{θ₁} sin θ dθ = R²·(π/2)·(cos θ₂ − cos θ₁)
area_ana = R_sphere^2 * (π/2) * (cos(θ₂) - cos(θ₁))
@printf("  Numerical (2×2 Gauss): %.6f\n", area_num)
@printf("  Analytical:            %.6f\n", area_ana)
@printf("  Relative error:        %.2e\n",  abs(area_num - area_ana)/area_ana)
