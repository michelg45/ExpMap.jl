```@meta
CurrentModule = ExpMap
```

# T\_SE3\_input\_data — Full data for the SE(3) tangent operator

`T_SE3_input_data` computes **all** quantities stored in [`T_SE3_data`](@ref)
in a single shared evaluation: first and second derivatives of `T(φ)` and
`T⁻¹(φ)` with respect to `φ`, together with the derived translation vectors.
This is the entry point to use when both the operator and its Gâteaux
derivative are needed.

## T\_SE3\_input\_data background

Given `p = [u; φ]`, the computation proceeds in layers:

```
Layer 0 — scalar coefficients
   x = ‖φ‖,  φ̂ = φ/x
   α, α′, α″, β, β′, β″  from T_functions(x; der=2)
   γ, γ′                  from invT_functions(x; der=1)

Layer 1 — W and its derivatives
   W   = φ̂⊗φ̂ − I,          dφ̂ = (I − φ̂⊗φ̂)/x
   dW[i] = ∂W/∂φᵢ

Layer 2 — first derivatives of T and T⁻¹
   dT[i]    = ∂T/∂φᵢ
   dinvT[i] = ∂T⁻¹/∂φᵢ

Layer 3 — second derivatives (lower triangle, then symmetrised)
   d²φ̂[k]  = ∂²φ̂/∂φₖ∂φⱼ
   d²W[i,j] = ∂²W/∂φᵢ∂φⱼ
   d²T[i,j] = ∂²T/∂φᵢ∂φⱼ

Layer 4 — translation contributions
   dTu[:,j]    = dT[j]·u     (assembled as MAT3 columns)
   d²Tu[i][:,j]  = d²T[i,j]·u
   d²Ttu[i][:,j] = (d²T[i,j])ᵀ·u
```

## T\_SE3\_input\_data small-angle limit

For `‖φ‖² < 1e-8`:

```
T ≈ I − ½φ̃,   T⁻¹ ≈ I + ½φ̃
dT[i] = −½ tilde(eᵢ),   dinvT[i] = +½ tilde(eᵢ)
d²Tu[i] = d²Ttu[i] = 0₃ₓ₃
```

**Dependencies:** [`T_functions`](@ref), [`invT_functions`](@ref),
[`tilde`](@ref), [`outerp`](@ref)

**See also:** [`T_SE3_data`](@ref), [`T_SE3`](@ref), [`DT_SE3`](@ref),
[`DinvT_SE3`](@ref)

```@docs
T_SE3_input_data
```
