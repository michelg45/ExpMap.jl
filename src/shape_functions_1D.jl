# ==============================================================================
# shape_functions_1D.jl — 1D Lagrange shape functions and their derivatives
# ==============================================================================
#
# Evaluates the Lagrange interpolation polynomials and their derivatives at a
# single point ξ on the reference segment [−1, +1].
#
# CALL
#   F, DF = shape_functions_1D(Nnodes, ξ)
#
# INPUT
#   Nnodes :: Int     number of nodes: 2, 3 or 4
#   ξ      :: Float64 evaluation point on [−1, +1]
#
# OUTPUT
#   F  :: Vector{Float64}   shape function values at ξ,   length = Nnodes
#   DF :: Vector{Float64}   dF/dξ values at ξ,            length = Nnodes
#
# ELEMENTS AND NODE LOCATIONS
#
#   Nnodes = 2   linear      L2   nodes at  ξ = −1, +1
#   Nnodes = 3   quadratic   L3   nodes at  ξ = −1,  0, +1
#   Nnodes = 4   cubic       L4   nodes at  ξ = −1, −1/3, +1/3, +1
#
# PARTITION OF UNITY
#   ∑ F[i] = 1  and  ∑ DF[i] = 0  for all ξ.
#
# ==============================================================================

"""
    shape_functions_1D(Nnodes::Int, ξ::Float64) → (F, DF)

Evaluate the 1D Lagrange shape functions and their ξ-derivatives at point `ξ`
on the reference segment [−1, +1].

- `Nnodes = 2`: linear element L2, nodes at ξ = −1, +1.
- `Nnodes = 3`: quadratic element L3, nodes at ξ = −1, 0, +1.
- `Nnodes = 4`: cubic element L4, nodes at ξ = −1, −1/3, +1/3, +1.

Returns `F::Vector{Float64}` (length `Nnodes`) and `DF::Vector{Float64}` (length `Nnodes`).
"""
function shape_functions_1D(Nnodes::Int, ξ::Float64)

    if Nnodes == 2
        # Linear L2: N₁ = ½(1−ξ),  N₂ = ½(1+ξ)
        F  = [0.5*(1.0 - ξ);  0.5*(1.0 + ξ)]
        DF = [-0.5;  0.5]

    elseif Nnodes == 3
        # Quadratic L3: nodes at −1, 0, +1
        F  = [0.5*ξ*(ξ - 1.0);  1.0 - ξ*ξ;  0.5*ξ*(1.0 + ξ)]
        DF = [ξ - 0.5;  -2.0*ξ;  ξ + 0.5]

    elseif Nnodes == 4
        # Cubic L4: nodes at −1, −1/3, +1/3, +1
        # Factor 9/16 = 0.5625 arises from the Lagrange product formula
        F  = 0.5625 * [-(1.0 - ξ)*(1.0/3.0 + ξ)*(1.0/3.0 - ξ);
                        3.0*(1.0 - ξ)*(1.0 + ξ)*(1.0/3.0 - ξ);
                        3.0*(1.0 - ξ)*(1.0 + ξ)*(1.0/3.0 + ξ);
                       -(1.0 + ξ)*(1.0/3.0 + ξ)*(1.0/3.0 - ξ)]
        DF = [-1.6875*ξ*ξ + 1.125*ξ + 0.0625;
               5.0625*ξ*ξ - 1.125*ξ - 1.6875;
              -5.0625*ξ*ξ - 1.125*ξ + 1.6875;
               1.6875*ξ*ξ + 1.125*ξ - 0.0625]
    end

    return F, DF
end
