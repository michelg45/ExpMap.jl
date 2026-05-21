# ==============================================================================
# sphere_Q2s_plot.jl — GLMakie visualisation of a spherical Q2s element
# ==============================================================================
#
# Sphere  R = 10, centre C = [2, 0, 0].
# Patch   z ∈ [6, 8],  φ ∈ [0, π/2].
# 8-node serendipity element (Q2s) defined by 8 NodeFrames.
#
# The script:
#   1. builds the 8 nodal frames via SE(3) interpolation
#   2. evaluates frame_interpol_2D on a 30×30 parametric grid
#   3. displays the z-coloured surface, the parametric mesh lines, and
#      the local frame axes at each node
# ==============================================================================

if !@isdefined(ExpMap)
    include(joinpath(@__DIR__, "..", "src", "ExpMap.jl"))
    using .ExpMap
end
using GLMakie

# ── Sphere parameters ─────────────────────────────────────────────────────────
const R_sphere = 10.0
const C_sphere = VEC3(2.0, 0.0, 0.0)

# ── Local frame on the sphere at (θ, φ) ───────────────────────────────────────
function sphere_frame(θ::Float64, φ::Float64)
    sθ, cθ = sin(θ), cos(θ)
    sφ, cφ = sin(φ), cos(φ)
    x   = C_sphere + R_sphere * VEC3(sθ*cφ, sθ*sφ, cθ)
    e_x = VEC3(-sφ,      cφ,     0.0)   # horizontal tangent  e_φ
    e_y = VEC3(-cθ*cφ, -cθ*sφ,  sθ )   # meridional tangent −e_θ
    e_z = VEC3( sθ*cφ,  sθ*sφ,  cθ )   # outward normal      n̂
    return NodeFrame(x, invR_SO3(MAT3(e_x, e_y, e_z)))
end

# ── 8 nodal frames ────────────────────────────────────────────────────────────
θ₁  = acos(3.0/5.0)   # latitude circle z = 6
θ₂  = acos(4.0/5.0)   # latitude circle z = 8
θ_m = π/4             # mid-latitude angle  (θ₁ + θ₂ = π/2)

H = Vector{NodeFrame}(undef, 8)
H[1] = sphere_frame(θ₁,  0.0 )
H[2] = sphere_frame(θ₁,  π/2 )
H[3] = sphere_frame(θ₂,  π/2 )
H[4] = sphere_frame(θ₂,  0.0 )
H[5] = sphere_frame(θ₁,  π/4 )
H[6] = sphere_frame(θ_m, π/2 )
H[7] = sphere_frame(θ₂,  π/4 )
H[8] = sphere_frame(θ_m, 0.0 )

# ── 30×30 mesh via SE(3) interpolation ───────────────────────────────────────
println("Building 30×30 mesh...")
N   = 30
ξv  = LinRange(-1.0, 1.0, N+1)   # 31 values of ξ₁ and ξ₂

X  = zeros(Float32, N+1, N+1)
Y  = zeros(Float32, N+1, N+1)
Z  = zeros(Float32, N+1, N+1)
Nx = zeros(Float32, N+1, N+1)
Ny = zeros(Float32, N+1, N+1)
Nz = zeros(Float32, N+1, N+1)

for j = 1:N+1, i = 1:N+1
    ξ = [ξv[i], ξv[j]]
    F, DF = shape_functions_2D(8, ξ)
    H_P, _, g_P, n_P, _ = frame_interpol_2D(H, F, DF)
    X[i,j]  = H_P.x[1];  Nx[i,j] = n_P[1]
    Y[i,j]  = H_P.x[2];  Ny[i,j] = n_P[2]
    Z[i,j]  = H_P.x[3];  Nz[i,j] = n_P[3]
end
println("  done — $(N*N) quadrilaterals, $((N+1)^2) nodes.")

# ── GLMakie figure ────────────────────────────────────────────────────────────
fig = Figure(size = (1100, 800))

ax = Axis3(fig[1, 1];
    title   = "Spherical Q2s element — SE(3) interpolation, $(N)×$(N) mesh",
    xlabel  = "x", ylabel = "y", zlabel = "z",
    aspect  = :data,
    azimuth = π/4,
    elevation = π/8)

# ── Coloured surface (colour = z altitude) ────────────────────────────────────
sf = surface!(ax, X, Y, Z;
    color    = Z,
    colormap = :viridis,
    shading  = true)

Colorbar(fig[1, 2], sf; label = "z", width = 16, tellheight = false)

# ── Parametric mesh lines ─────────────────────────────────────────────────────
# ξ₁ = const (φ-curves)
for i = 1:N+1
    lines!(ax, X[i, :], Y[i, :], Z[i, :];
        color = (:white, 0.55), linewidth = 0.7)
end
# ξ₂ = const (θ-curves)
for j = 1:N+1
    lines!(ax, X[:, j], Y[:, j], Z[:, j];
        color = (:white, 0.55), linewidth = 0.7)
end

# ── Local frames at the 8 nodes ───────────────────────────────────────────────
s      = 0.9f0                                    # arrow length
clrs   = [:dodgerblue, :limegreen, :orangered2]   # e_x, e_y, e_z

for k = 1:8
    xk = H[k].x
    R  = R_SO3(H[k].phi)
    p  = Point3f(xk[1], xk[2], xk[3])

    # Node marker
    scatter!(ax, [p]; color = :white, markersize = 12, strokecolor = :black,
             strokewidth = 1)

    # Three local axes
    for col = 1:3
        d = Vec3f(R[1, col], R[2, col], R[3, col])
        arrows3d!(ax, [p], [s * d];
            color       = clrs[col],
            shaftradius = 0.04,
            tipradius   = 0.10,
            tiplength   = 0.22)
    end
end

# ── Legend ────────────────────────────────────────────────────────────────────
leg_elems = [LineElement(color = c, linewidth = 3) for c in clrs]
leg_labels = ["e_x = eφ  (horizontal tangent)",
              "e_y = −eθ (meridional tangent)",
              "e_z = n̂   (outward normal)"]
Legend(fig[2, 1], leg_elems, leg_labels, "Local frame axes";
    orientation = :horizontal, tellwidth = false,
    framevisible = false)

display(fig)
println("GLMakie window opened.")
save(joinpath(@__DIR__,"shell.png"),fig)
