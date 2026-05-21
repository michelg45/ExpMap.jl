```@meta
CurrentModule = ExpMap
```

# frame\_interpol\_1D — SE(3) Frame Interpolation on a 1D Element

`frame_interpol_1D` interpolates a frame curve from a set of `Nnodes` nodal
frames `H[k] ∈ SE(3)` at a given point `ξ`, using Lagrange shape functions.
It returns the interpolated frame, the generalised strain vector, and the
nodal sensitivity matrices needed for tangent stiffness assembly.

## frame\_interpol\_1D algorithm

The interpolated frame `H_P` is the **SE(3) weighted mean** of the nodal
frames, defined as the fixed point of the iterative equation:

```
H_P = argmin  ½ Σ F[k] ‖log_SE3(H_P⁻¹·H[k])‖²
```

It is found by the Newton iteration:

```
p[k]  ←  log_SE3(H_P⁻¹·H[k])         relative motion in se(3)
r     ←  Σ F[k]·p[k]                  weighted residual
A     ←  Σ F[k]·invT_SE3(−p[k])       tangent accumulation matrix
H_P   ←  H_P · exp_SE3(A⁻¹·r)         SE(3) update step
```

Convergence is declared when `‖r‖ < 1e-12 · ‖r₀‖` (relative tolerance
normalised by the first residual).

## frame\_interpol\_1D strain vector

Once `H_P` is converged, the generalised strain vector is assembled from the
scaled shape function derivatives `DFs[k] = (2/ℓ)·DF[k]`:

```
f_P = A⁻¹ · Σ DFs[k]·p[k]
```

The factor `2/ℓ` maps derivatives from the reference segment [−1,+1] to
physical arc length `s ∈ [0,ℓ]`.  Components 1–3 of `f_P` are
displacement-gradient measures; components 4–6 are curvatures.

## frame\_interpol\_1D sensitivity matrices

The nodal sensitivity `QN[k] = ∂f_P/∂p_k` is computed by the chain rule.
For each node `k`:

```
invTm  =  invT_SE3(a_k; trp=true)           (T⁻¹)ᵀ-based operator
invTp  =  invT_SE3(a_k; trp=false)          standard T⁻¹-based operator
DTmpf  =  DinvT_SE3(a_k, f_P; sign_p="-")  Gâteaux derivative
Q[k]   =  F[k]·invTm
B      =  F[k]·DTmpf + DFs[k]·I₆
QN[k]  =  B·invTp
C     -=  B·invTm
```

followed by the left-multiplication `A⁻¹` pass:

```
Q[k]   ←  A⁻¹·Q[k]
QN[k]  ←  A⁻¹·(QN[k] + C·Q[k])
```

## frame\_interpol\_1D example

```julia
# Two-node linear beam element
H1 = NodeFrame(VEC3(0.0, 0.0, 0.0), RV3(0.0, 0.0, 0.0))
H2 = NodeFrame(VEC3(1.0, 0.0, 0.0), RV3(0.0, 0.2, 0.0))
H  = [H1, H2]
ℓ  = 1.0

# Evaluate at the midpoint ξ = 0
F,  DF  = shape_functions_1D(2, 0.0)
H_P, f_P, QN = frame_interpol_1D(H, F, DF, ℓ)
```

**Dependencies:** [`log_SE3`](@ref), [`exp_SE3`](@ref), [`invT_SE3`](@ref),
[`T_SE3_input_data`](@ref), [`DinvT_SE3`](@ref), [`shape_functions_1D`](@ref)

**See also:** [`NodeFrame`](@ref), [`T_SE3_data`](@ref)

```@docs
frame_interpol_1D
```
