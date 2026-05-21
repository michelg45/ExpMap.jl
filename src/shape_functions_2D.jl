# ==============================================================================
# shape_functions_2D.jl — 2D Lagrange shape functions and their derivatives
# ==============================================================================
#
# Evaluates the Lagrange interpolation polynomials and their partial derivatives
# at a single point ξ = [ξ₁, ξ₂] on a 2D reference domain.
#
# CALL
#   F, DF = shape_functions_2D(Nnodes, ξ)
#
# INPUT
#   Nnodes :: Int              number of nodes: 3, 4, 5, 6, 8 or 9
#   ξ      :: Vector{Float64}  evaluation point [ξ₁, ξ₂]
#
# OUTPUT
#   F  :: Vector{Float64}     shape function values at ξ,  length = Nnodes
#   DF :: Matrix{Float64}     partial derivatives at ξ,    size   = Nnodes × 2
#                             DF[i,1] = ∂Fᵢ/∂ξ₁,  DF[i,2] = ∂Fᵢ/∂ξ₂
#
# SUPPORTED ELEMENTS
#
#  Reference triangle   vertices (0,0), (1,0), (0,1)
#                       natural coordinates ξ₁ ≥ 0, ξ₂ ≥ 0, ξ₁+ξ₂ ≤ 1
#
#   Nnodes = 3   linear triangle T3
#                nodes:  1=(0,0)  2=(1,0)  3=(0,1)
#
#   Nnodes = 6   quadratic triangle T6 (Lagrange)
#                nodes:  1=(0,0)  2=(1,0)  3=(0,1)
#                        4=(½,0)  5=(½,½)  6=(0,½)
#
#  Reference square   [−1,+1]²,  nodes labelled counter-clockwise
#
#   Nnodes = 4   bilinear quadrilateral Q1
#                corners: 1=(−1,−1)  2=(+1,−1)  3=(+1,+1)  4=(−1,+1)
#
#   Nnodes = 5   enriched Q1 with central bubble node
#                corners 1–4 as Q1; node 5 = (0,0) (central node)
#                Note: the bubble function f_E = (1−ξ₁)²(1−ξ₂)² vanishes
#                on edges ξ₁=+1 and ξ₂=+1 but not on ξ₁=−1, ξ₂=−1.
#
#   Nnodes = 8   serendipity quadrilateral Q2s (8-node, no centre node)
#                corners 1–4 as Q1; mid-side nodes:
#                5=(0,−1)  6=(+1,0)  7=(0,+1)  8=(−1,0)
#
#   Nnodes = 9   biquadratic quadrilateral Q2 (Lagrange, 9-node)
#                corners 1–4 as Q1; mid-side nodes 5–8 as Q2s;
#                centre node 9=(0,0)
#
# PARTITION OF UNITY
#   ∑ F[i] = 1  and  ∑ DF[i,:] = [0 0]  for all ξ.
#
# ==============================================================================

"""
    shape_functions_2D(Nnodes::Int, ξ::Vector{Float64}) → (F, DF)

Evaluate 2D Lagrange shape functions and their partial derivatives at point `ξ = [ξ₁, ξ₂]`.

Supported elements:

| `Nnodes` | Element | Reference domain |
|----------|---------|-----------------|
| 3  | linear triangle T3            | triangle (0,0)–(1,0)–(0,1) |
| 4  | bilinear quad Q1              | square [−1,+1]² |
| 5  | enriched Q1 + bubble node     | square [−1,+1]² |
| 6  | quadratic triangle T6         | triangle (0,0)–(1,0)–(0,1) |
| 8  | serendipity quad Q2s          | square [−1,+1]² |
| 9  | biquadratic quad Q2 (Lagrange)| square [−1,+1]² |

Returns `F::Vector{Float64}` (length `Nnodes`) and `DF::Matrix{Float64}` (`Nnodes×2`),
where `DF[i,j] = ∂Fᵢ/∂ξⱼ`.
"""
function shape_functions_2D(Nnodes::Int, ξ::Vector{Float64})

    if Nnodes == 3
        # Linear triangle T3 — natural coordinates on (0,0),(1,0),(0,1)
        F  = [1.0 - ξ[1] - ξ[2];  ξ[1];  ξ[2]]

        DF = [-1.0  -1.0;
               1.0   0.0;
               0.0   1.0]

    elseif Nnodes == 4
        # Bilinear quadrilateral Q1 — square [−1,+1]²
        F  = [0.25*(1.0 - ξ[1])*(1.0 - ξ[2]);
              0.25*(1.0 + ξ[1])*(1.0 - ξ[2]);
              0.25*(1.0 + ξ[1])*(1.0 + ξ[2]);
              0.25*(1.0 - ξ[1])*(1.0 + ξ[2])]

        DF = [-0.25*(1.0 - ξ[2])   -0.25*(1.0 - ξ[1]);
               0.25*(1.0 - ξ[2])   -0.25*(1.0 + ξ[1]);
               0.25*(1.0 + ξ[2])    0.25*(1.0 + ξ[1]);
              -0.25*(1.0 + ξ[2])    0.25*(1.0 - ξ[1])]

    elseif Nnodes == 5
        # Enriched Q1 with central bubble node — square [−1,+1]²
        # Bubble function f_E = (1−ξ₁)²(1−ξ₂)² vanishes on edges ξ₁=+1 and ξ₂=+1
        f_E    = (1.0 - ξ[1])^2 * (1.0 - ξ[2])^2
        d1_f_E = -2.0*ξ[1] * (1.0 - ξ[2])^2
        d2_f_E = -2.0*ξ[2] * (1.0 - ξ[1])^2

        F  = [0.25*(1.0 - ξ[1])*(1.0 - ξ[2]) - 0.25*f_E;
              0.25*(1.0 + ξ[1])*(1.0 - ξ[2]) - 0.25*f_E;
              0.25*(1.0 + ξ[1])*(1.0 + ξ[2]) - 0.25*f_E;
              0.25*(1.0 - ξ[1])*(1.0 + ξ[2]) - 0.25*f_E;
              f_E]

        DF = [-0.25*(1.0 - ξ[2]) - 0.25*d1_f_E   -0.25*(1.0 - ξ[1]) - 0.25*d2_f_E;
               0.25*(1.0 - ξ[2]) - 0.25*d1_f_E   -0.25*(1.0 + ξ[1]) - 0.25*d2_f_E;
               0.25*(1.0 + ξ[2]) - 0.25*d1_f_E    0.25*(1.0 + ξ[1]) - 0.25*d2_f_E;
              -0.25*(1.0 + ξ[2]) - 0.25*d1_f_E    0.25*(1.0 - ξ[1]) - 0.25*d2_f_E;
               d1_f_E                              d2_f_E]

    elseif Nnodes == 6
        # Quadratic triangle T6 — natural coordinates on (0,0),(1,0),(0,1)
        # Mid-side nodes at (½,0), (½,½), (0,½)
        F  = [1.0 - 3.0*ξ[1] - 3.0*ξ[2] + 2.0*ξ[1]^2 + 4.0*ξ[1]*ξ[2] + 2.0*ξ[2]^2;
             -ξ[1] + 2.0*ξ[1]^2;
             -ξ[2] + 2.0*ξ[2]^2;
              4.0*ξ[1] - 4.0*ξ[1]^2 - 4.0*ξ[1]*ξ[2];
              4.0*ξ[1]*ξ[2];
              4.0*ξ[2] - 4.0*ξ[2]^2 - 4.0*ξ[1]*ξ[2]]

        DF = [-3.0 + 4.0*ξ[1] + 4.0*ξ[2]   -3.0 + 4.0*ξ[1] + 4.0*ξ[2];
              -1.0 + 4.0*ξ[1]                0.0;
               0.0                           -1.0 + 4.0*ξ[2];
               4.0 - 8.0*ξ[1] - 4.0*ξ[2]   -4.0*ξ[1];
               4.0*ξ[2]                       4.0*ξ[1];
              -4.0*ξ[2]                       4.0 - 8.0*ξ[2] - 4.0*ξ[1]]

    elseif Nnodes == 8
        # Serendipity quadrilateral Q2s — square [−1,+1]²
        # Corners 1–4 counter-clockwise; mid-side nodes 5–8
        η = ξ[1];  μ = ξ[2]

        F  = [-0.25*(1.0 - η)*(1.0 - μ)*(1.0 + η + μ);
              -0.25*(1.0 + η)*(1.0 - μ)*(1.0 - η + μ);
              -0.25*(1.0 + η)*(1.0 + μ)*(1.0 - η - μ);
              -0.25*(1.0 - η)*(1.0 + μ)*(1.0 + η - μ);
               0.5*(1.0 - η^2)*(1.0 - μ);
               0.5*(1.0 + η)*(1.0 - μ^2);
               0.5*(1.0 - η^2)*(1.0 + μ);
               0.5*(1.0 - η)*(1.0 - μ^2)]

        DF = [-0.25*(1.0 - μ)*(-μ - 2.0*η)   -0.25*(1.0 - η)*(-η - 2.0*μ);
              -0.25*(1.0 - μ)*( μ - 2.0*η)   -0.25*(1.0 + η)*( η - 2.0*μ);
              -0.25*(1.0 + μ)*(-μ - 2.0*η)   -0.25*(1.0 + η)*(-η - 2.0*μ);
              -0.25*(1.0 + μ)*( μ - 2.0*η)   -0.25*(1.0 - η)*( η - 2.0*μ);
              -(1.0 - μ)*η                    -0.5*(1.0 - η^2);
               0.5*(1.0 - μ^2)               -(1.0 + η)*μ;
              -(1.0 + μ)*η                     0.5*(1.0 - η^2);
              -0.5*(1.0 - μ^2)               -(1.0 - η)*μ]

    elseif Nnodes == 9
        # Biquadratic quadrilateral Q2 (Lagrange, 9-node) — square [−1,+1]²
        # Corners 1–4, mid-side nodes 5–8, centre node 9
        η = ξ[1];  μ = ξ[2]

        # Auxiliary factors (reused across F and DF)
        Omη  = 1.0 - η;   Opη  = 1.0 + η
        Omμ  = 1.0 - μ;   Opμ  = 1.0 + μ
        Omη2 = 1.0 - η^2; Omμ2 = 1.0 - μ^2
        Om2η = 1.0 - 2.0*η; Op2η = 1.0 + 2.0*η
        Om2μ = 1.0 - 2.0*μ; Op2μ = 1.0 + 2.0*μ

        F  = [ 0.25*Omη*Omμ*η*μ;
              -0.25*Opη*Omμ*η*μ;
               0.25*Opη*Opμ*η*μ;
              -0.25*Omη*Opμ*η*μ;
              -0.5*μ*Omμ*Omη2;
               0.5*η*Opη*Omμ2;
               0.5*μ*Opμ*Omη2;
              -0.5*η*Omη*Omμ2;
               Omη2*Omμ2]

        DF = [ 0.25*μ*Omμ*Om2η   0.25*Omη*η*Om2μ;
              -0.25*μ*Omμ*Op2η  -0.25*Opη*η*Om2μ;
               0.25*μ*Opμ*Op2η   0.25*Opη*η*Op2μ;
              -0.25*μ*Opμ*Om2η  -0.25*Omη*η*Op2μ;
               μ*Omμ*η          -0.5*Omη2*Om2μ;
               0.5*Omμ2*Op2η    -η*Opη*μ;
              -μ*Opμ*η           0.5*Omη2*Op2μ;
              -0.5*Omμ2*Om2η     η*Omη*μ;
              -2.0*η*Omμ2       -2.0*μ*Omη2]
    end

    return F, DF
end
