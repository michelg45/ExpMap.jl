```@meta
CurrentModule = ExpMap
```

# frame\_interpol\_2D — SE(3) Frame Interpolation on a 2D Surface Element

`frame_interpol_2D` provides two methods for interpolating frame curves on 2D
surface elements using the SE(3) weighted mean.  Both share the same Newton
iteration as [`frame_interpol_1D`](@ref).

## frame\_interpol\_2D — Method 1: reference surface

Called without a precomputed Jacobian, this method interpolates the **reference
(undeformed) surface** and returns the geometric quantities needed to set up
the integration on that surface.

```
frame_interpol_2D(H, F, DF) → (H_P, f_P, g_P, n_P, J_P)
```

**Strain assembly (reference coordinates):**
```
f_P_ref[α] = invA⁻¹ · Σ_k DF[k,α]·p[k]     α = 1, 2
```

**Surface Jacobian** (2×2 matrix mapping reference ξ to physical in-plane coords):
```
J_P[i,α] = f_P_ref[α][i]     i, α = 1, 2
```

**Physical strain vectors** (after removing Jacobian distortion):
```
f_P[α] = Σ_β  f_P_ref[β] · J_P⁻¹[β,α]
```

**Tangent vectors and unit normal:**
```
g_P[α] = R(H_P.phi) · VEC3(f_P[α][1:3])
n_P    = (g_P[1] × g_P[2]) / (‖g_P[1]‖·‖g_P[2]‖)
```

The returned `J_P` is passed to Method 2 for all subsequent quadrature points
in the current configuration.

## frame\_interpol\_2D — Method 2: current surface with sensitivity

Called with the reference Jacobian `J` from Method 1, this method interpolates
the **current (instantaneous) surface** and computes nodal sensitivities.

```
frame_interpol_2D(H, F, DF, J) → (H_P, f_P, g_P, n_P, QN)
```

**Derivative scaling** (reference → physical coordinates):
```
DFs = DF · J⁻¹     (Nnodes × 2 matrix)
```

**Physical strain assembly:**
```
f_P[α] = invA⁻¹ · Σ_k DFs[k,α]·p[k]
```

**Sensitivities** `QN[α,k] = ∂f_P[α]/∂p_k` (6×6 matrix, `α ∈ {1,2}`):
Same chain-rule structure as [`frame_interpol_1D`](@ref), applied independently
to each surface direction:
```
invTm    = invT_SE3(a_k; trp=true)
invTp    = invT_SE3(a_k; trp=false)
DTmpf[α] = DinvT_SE3(a_k, f_P[α]; sign_p="-")
Q[k]     = F[k]·invTm
B[α]     = F[k]·DTmpf[α] + DFs[k,α]·I₆
QN[α,k]  = B[α]·invTp
C[α]    -= B[α]·invTm
```
followed by:
```
Q[k]    ← invA⁻¹·Q[k]
QN[α,k] ← invA⁻¹·(QN[α,k] + C[α]·Q[k])
```

!!! note "Verification status"
    Method 2 is provided as-is and has not yet been fully validated against an
    independent reference.  The sensitivity structure follows the same pattern
    as [`frame_interpol_1D`](@ref).

## frame\_interpol\_2D example

```julia
# Reference surface setup (once per element)
H = [NodeFrame(...) for k = 1:Nnodes]
ξ_ref, w_ref = gauss_points(2, 4)

F, DF = shape_functions_2D(Nnodes, ξ_ref[1])
H_P_ref, f_P_ref, g_P_ref, n_P_ref, J_P = frame_interpol_2D(H, F, DF)

# Current configuration at the same point
F, DF = shape_functions_2D(Nnodes, ξ_ref[1])
H_P, f_P, g_P, n_P, QN = frame_interpol_2D(H_current, F, DF, J_P)
```

**Dependencies:** [`log_SE3`](@ref), [`exp_SE3`](@ref), [`invT_SE3`](@ref),
[`R_SO3`](@ref), [`T_SE3_input_data`](@ref), [`DinvT_SE3`](@ref),
[`shape_functions_2D`](@ref), [`gauss_points`](@ref)

**See also:** [`frame_interpol_1D`](@ref), [`NodeFrame`](@ref)

```@docs
frame_interpol_2D
```
