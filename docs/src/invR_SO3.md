```@meta
CurrentModule = ExpMap
```

# invR_SO3 — Logarithmic Map SO(3)

`invR_SO3` implements the **logarithmic map** from SO(3) back to its Lie
algebra so(3): given a 3×3 rotation matrix R ∈ SO(3), it recovers the
Cartesian Rotation Vector (CRV) ψ such that `R_SO3(ψ) = R`.

## Algorithm — Cayley / quaternion method

Any rotation matrix R can be decomposed via its quaternion representation.
The algorithm proceeds as follows:

1. Compute symmetric part **S = R + Rᵀ** and skew part **G = R − Rᵀ**.
2. Extract the quaternion vector part **e** component-wise
   (magnitude from S, sign from G).
3. Recover `cos(θ/2)` and `sin(θ/2) = ‖e‖`, then **θ = 2 atan(‖e‖, cos(θ/2))**.
4. Return `ψ = (θ / ‖e‖) · e`, or the first-order approximation `ψ ≈ 2e`
   for small θ.

!!! tip "Numerical superiority"
    This Cayley approach is numerically superior to the classical arccos/arcsin
    formula because it remains well-conditioned for **both small and large
    angles**, including θ near π.

## Precision

The small-angle threshold is `PREC = √eps(Float64) ≈ 1.49e-8` rad.
Below this threshold `sin(θ/2) ≈ θ/2` to machine precision, and the
first-order approximation `ψ ≈ 2e` is used.

## invR\_SO3 uniqueness

The result is unique for θ ∈ [0, π).  At θ = π the axis is determined only
up to sign; `invR_SO3` returns one valid choice.

## invR\_SO3 methods

### `invR_SO3(R::MAT3) → RV3`

Main logarithmic map function.

### `RV3(R::MAT3) → RV3`

Convenience alias allowing constructor-style syntax: `RV3(R)`.

**See also:** [`R_SO3`](@ref), [`rv3_comp_rule`](@ref)

```@docs
invR_SO3
```
