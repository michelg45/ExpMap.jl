# ==============================================================================
# curve_1D_plot.jl — GLMakie visualisation of SE(3) frame interpolation on 1D
# ==============================================================================
#
# Reads the 4 reference NodeFrames from curve_nodes.jl (L4 element on [-1,+1]),
# interpolates the frame curve at 101 uniformly-spaced parameter values, and
# opens two separate windows:
#
#   Figure 1 (single panel)
#     Curve + local frame axes (e_x, e_y, e_z) at the 4 reference nodes
#
#   Figure 2 (two panels side by side)
#     Left  — curve + displacement-gradient vectors  f_P[1:3]  at sample points
#     Right — curve + curvature vectors              f_P[4:6]  at sample points
#
# The strain vector f_P is the output of frame_interpol_1D scaled by 2/ℓ:
#   f_P[1:3]  translational strain (displacement gradient along the curve)
#   f_P[4:6]  rotational strain    (curvature / twist)
# ==============================================================================

if !@isdefined(ExpMap)
    include(joinpath(@__DIR__, "..", "src", "ExpMap.jl"))
    using .ExpMap
end
using GLMakie
using LinearAlgebra: norm
using Printf

# ── Load reference nodes ───────────────────────────────────────────────────────
# curve_nodes.jl defines:  Nnodes = 4,  H = Vector{NodeFrame}(undef, 4)
include(joinpath(@__DIR__, "curve_nodes.jl"))

# ── Element arc length (piecewise chord sum between consecutive nodes) ─────────
ℓ = sum(norm2(H[k+1].x - H[k].x) for k = 1:Nnodes-1)
@printf("Element arc length  ℓ ≈ %.4f\n", ℓ)

# ── Evaluate the interpolated curve at Nmesh+1 points ─────────────────────────
Nmesh = 100
ξv    = LinRange(-1.0, 1.0, Nmesh + 1)

pts      = Vector{Point3f}(undef, Nmesh + 1)   # interpolated positions
grad_vec = Vector{Vec3f}(undef, Nmesh + 1)     # f_P[1:3]  displacement gradient
curv_vec = Vector{Vec3f}(undef, Nmesh + 1)     # f_P[4:6]  curvature

for i = 1:Nmesh + 1
    F, DF       = shape_functions_1D(Nnodes, Float64(ξv[i]))
    H_P, f_P, _ = frame_interpol_1D(H, F, DF, ℓ)
    pts[i]      = Point3f(H_P.x[1],  H_P.x[2],  H_P.x[3])
    grad_vec[i] = Vec3f(f_P[1], f_P[2], f_P[3])
    curv_vec[i] = Vec3f(f_P[4], f_P[5], f_P[6])
end

# ── Adaptive arrow scaling (arrows ≈ 8% of element length) ────────────────────
max_grad = maximum(norm(v) for v in grad_vec)
max_curv = maximum(norm(v) for v in curv_vec)
s_node   = Float32(0.08 * ℓ)                        # frame axes at nodes
s_grad   = Float32(0.08 * ℓ / max(max_grad, 1e-6))  # displacement gradient
s_curv   = Float32(0.08 * ℓ / max(max_curv, 1e-6))  # curvature
@printf("  max |f_P[1:3]| = %.4f   arrow scale s_grad = %.4f\n", max_grad, s_grad)
@printf("  max |f_P[4:6]| = %.4f   arrow scale s_curv = %.4f\n", max_curv, s_curv)

# Sample every Nstep-th mesh point for arrow plots
Nstep     = 10
arrow_idx = 1:Nstep:Nmesh+1

# Shared arrow geometry (shaft and tip dimensions)
shaft_r = 0.025f0
tip_r   = 0.055f0
tip_len = 0.090f0

# Shared Axis3 keyword arguments
axis_kw = (
    xlabel         = "x", ylabel = "y", zlabel = "z",
    aspect         = :data,
    azimuth        = 1.1,
    elevation      = 0.35,
    xgridvisible   = false, ygridvisible   = false, zgridvisible   = false,
    xypanelvisible = false, yzpanelvisible = false, xzpanelvisible = false,
)

# ==============================================================================
# Figure 1 — local frame axes at the 4 reference nodes
# ==============================================================================
fig1 = Figure(size = (700, 620))

ax1 = Axis3(fig1[1, 1]; title = "Local frame axes at reference nodes", axis_kw...)
lines!(ax1, pts; color = :steelblue, linewidth = 3)

clrs_frame  = [:dodgerblue, :limegreen, :orangered2]
label_frame = ["e_x", "e_y", "e_z"]

for k = 1:Nnodes
    R = R_SO3(H[k].phi)
    p = Point3f(H[k].x[1], H[k].x[2], H[k].x[3])
    scatter!(ax1, [p]; color = :white, markersize = 14,
             strokecolor = :black, strokewidth = 1)
    for (col, c) in enumerate(clrs_frame)
        d = Vec3f(R[1, col], R[2, col], R[3, col])
        arrows3d!(ax1, [p], [s_node * d];
            color = c, shaftradius = shaft_r, tipradius = tip_r, tiplength = tip_len)
    end
end

leg_elems = [LineElement(color = c, linewidth = 3) for c in clrs_frame]
Legend(fig1[2, 1], leg_elems, label_frame, "Frame axes";
    orientation = :horizontal, tellwidth = false, framevisible = false)

rowsize!(fig1.layout, 2, Relative(0.08))
display(fig1)
println("Figure 1 opened.")

# ==============================================================================
# Figure 2 — displacement gradient (left) and curvature (right)
# ==============================================================================
fig2 = Figure(size = (1200, 620))

ax2 = Axis3(fig2[1, 1]; title = "Displacement gradient  f_P[1:3]", axis_kw...)
ax3 = Axis3(fig2[1, 2]; title = "Curvature  f_P[4:6]",             axis_kw...)

for ax in (ax2, ax3)
    lines!(ax, pts; color = :steelblue, linewidth = 3)
end

# ── Left panel — displacement-gradient vectors ─────────────────────────────────
arr_p2 = pts[arrow_idx]
arr_d2 = grad_vec[arrow_idx]
mags2  = Float32[norm(v) for v in arr_d2]
clim2  = (0f0, Float32(max_grad))

arrows3d!(ax2, arr_p2, s_grad .* arr_d2;
    color = mags2, colormap = :plasma, colorrange = clim2,
    shaftradius = shaft_r, tipradius = tip_r, tiplength = tip_len)

Colorbar(fig2[2, 1]; colormap = :plasma, limits = clim2,
    label = "|f_P[1:3]|", vertical = false, tellwidth = false)

# ── Right panel — curvature vectors ───────────────────────────────────────────
arr_p3 = pts[arrow_idx]
arr_d3 = curv_vec[arrow_idx]
mags3  = Float32[norm(v) for v in arr_d3]
clim3  = (0f0, Float32(max_curv))

arrows3d!(ax3, arr_p3, s_curv .* arr_d3;
    color = mags3, colormap = :inferno, colorrange = clim3,
    shaftradius = shaft_r, tipradius = tip_r, tiplength = tip_len)

Colorbar(fig2[2, 2]; colormap = :inferno, limits = clim3,
    label = "|f_P[4:6]|", vertical = false, tellwidth = false)

rowsize!(fig2.layout, 2, Relative(0.08))
display(fig2)
println("Figure 2 opened.")
save(joinpath(@__DIR__,"curve_frames.png"),fig1)
save(joinpath(@__DIR__,"curve_gradients.png"),fig2)