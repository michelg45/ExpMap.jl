# ==============================================================================
# gauss_points.jl — Gauss–Legendre quadrature points and weights
# ==============================================================================
#
# Returns the abscissae and weights for Gauss–Legendre quadrature rules on
# standard reference domains in 1D and 2D.
#
# CALL
#   x, w = gauss_points(ndim, npoints)
#
# INPUT
#   ndim    :: Int   spatial dimension: 1 (line) or 2 (quad / triangle)
#   npoints :: Int   number of quadrature points (see table below)
#
# OUTPUT
#   x  — quadrature abscissae
#          ndim == 1  :  Vector{Float64} of length npoints
#          ndim == 2  :  Vector{Vector{Float64}} of length npoints,
#                        each inner vector is a 2-component coordinate
#   w  — quadrature weights  (Vector{Float64} of length npoints)
#
# REFERENCE DOMAINS AND SUPPORTED RULES
#
#   ndim == 1   segment [−1, +1]
#     npoints = 1   1-point rule   (exact for polynomials up to degree 1)
#     npoints = 2   2-point rule   (exact up to degree 3)
#     npoints = 3   3-point rule   (exact up to degree 5)
#     npoints = 4   4-point rule   (exact up to degree 7)
#
#   ndim == 2, quad   square [−1,+1]²  (tensor-product rules)
#     npoints = 1   1-point rule   (2×2 zeroth-order Gauss)
#     npoints = 4   4-point rule   (2×2 Gauss–Legendre)
#     npoints = 9   9-point rule   (3×3 Gauss–Legendre)
#
#   ndim == 2, triangle  reference triangle with vertices (0,0),(1,0),(0,1)
#     npoints = 7   7-point rule (Dunavant, degree 5)
#
# NOTE
#   The weights for the 1D rules sum to 2 (length of [−1,+1]).
#   The weights for the quad rules sum to 4 (area of [−1,+1]²).
#   The weights for the 7-point triangular rule sum to 1/2 (area of the
#   reference triangle), after the 1/2 factor applied to the Dunavant weights.
#
# ==============================================================================

"""
    gauss_points(ndim::Int, npoints::Int) → (x, w)

Return the abscissae `x` and weights `w` for a Gauss–Legendre quadrature rule
on the standard reference domain.

- `ndim = 1`: segment [−1, +1]; `npoints` ∈ {1, 2, 3, 4}.
  `x` and `w` are `Vector{Float64}`.

- `ndim = 2`, quad [−1,+1]²: `npoints` ∈ {1, 4, 9}.
  `x` is a `Vector{Vector{Float64}}`, `w` is a `Vector{Float64}`.

- `ndim = 2`, triangle with vertices (0,0),(1,0),(0,1): `npoints = 7`
  (Dunavant degree-5 rule). `x` is a `Vector{Vector{Float64}}`,
  `w` is a `Vector{Float64}` (weights scaled by 1/2, summing to area = 1/2).
"""
function gauss_points(ndim::Int, npoints::Int)

    # ── 1D Gauss–Legendre rules on [−1, +1] ──────────────────────────────────
    if ndim == 1

        if npoints == 1
            # Midpoint rule — exact for degree ≤ 1
            x = [0.0]
            w = [2.0]

        elseif npoints == 2
            # 2-point rule — exact for degree ≤ 3
            x = 1.0/sqrt(3.0) * [-1.0; 1.0]
            w = ones(2)

        elseif npoints == 3
            # 3-point rule — exact for degree ≤ 5
            x = sqrt(3.0/5.0) * [0.0; -1.0; 1.0]
            w = 1.0/9.0 * [8.0; 5.0; 5.0]

        elseif npoints == 4
            # 4-point rule — exact for degree ≤ 7
            a = sqrt(3.0/7.0 - 2.0/7.0*sqrt(6.0/5.0))
            b = sqrt(3.0/7.0 + 2.0/7.0*sqrt(6.0/5.0))
            x = [-a; a; -b; b]
            w = 1.0/36.0 * (18.0*ones(4) + sqrt(30.0)*[1.0; 1.0; -1.0; -1.0])
        end

    # ── 2D rules ─────────────────────────────────────────────────────────────
    elseif ndim == 2

        if npoints == 1
            # 1-point rule on [−1,+1]² — centroid
            x = [[0.0; 0.0]]
            w = [4.0]

        elseif npoints == 4
            # 2×2 Gauss–Legendre on [−1,+1]²
            c = 1.0/sqrt(3.0)
            x = [[-c, -c], [c, -c], [-c, c], [c, c]]
            w = ones(4)

        elseif npoints == 7
            # Dunavant degree-5 rule on the reference triangle
            # with vertices (0,0), (1,0), (0,1).
            # Barycentric coordinates (λ₁, λ₂, λ₃) are mapped to 2D via
            #   [x; y] = λ₂·[1;0] + λ₃·[0;1]
            # (λ₁ = 1 − λ₂ − λ₃ is the coordinate associated with vertex (0,0))
            # Weights are pre-scaled by 1/2 so that Σw = area(triangle) = 1/2.

            alpha_1 = 0.059715871789770;  beta_1 = 0.470142064105115
            alpha_2 = 0.797426985353087;  beta_2 = 0.101286507323456

            # Barycentric coordinates of the 7 quadrature points
            bary = Vector{Vector{Float64}}(undef, 7)
            bary[1] = [1.0/3.0, 1.0/3.0, 1.0/3.0]          # centroid
            bary[2] = [alpha_1, beta_1,  beta_1 ]
            bary[3] = [beta_1,  alpha_1, beta_1 ]
            bary[4] = [beta_1,  beta_1,  alpha_1]
            bary[5] = [alpha_2, beta_2,  beta_2 ]
            bary[6] = [beta_2,  alpha_2, beta_2 ]
            bary[7] = [beta_2,  beta_2,  alpha_2]

            # Vertex coordinates of the reference triangle (columns)
            #   v₁ = (0,0),  v₂ = (1,0),  v₃ = (0,1)
            vertices = zeros(2, 3)
            vertices[1, 2] = 1.0
            vertices[2, 3] = 1.0

            # Cartesian coordinates: [x;y] = vertices * λ
            x = [vertices * bary[i] for i = 1:7]

            # Dunavant weights (sum = 1 before scaling)
            w = zeros(7)
            w[1]   = 0.225
            w[2:4] .= 0.132394152788506
            w[5:7] .= 0.125939180544827
            w .*= 0.5      # scale so that Σw = area of reference triangle = 1/2

        elseif npoints == 9
            # 3×3 Gauss–Legendre on [−1,+1]²
            c = sqrt(3.0/5.0)
            x = [[-c; -c]; [[0.0; -c]]; [[c; -c]];
                 [[-c;  0.0]]; [zeros(2)]; [[c;  0.0]];
                 [[-c;  c]]; [[0.0;  c]]; [[c;  c]]]
            w = 1/81.0 * [25.0; 40.0; 25.0; 40.0; 64.0; 40.0; 25.0; 40.0; 25.0]
        end

    end

    return x, w
end
