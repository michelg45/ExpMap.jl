```@meta
CurrentModule = ExpMap
```

# exp_SE3 — Exponential Map SE(3)

`exp_SE3` implements the **exponential map** from the Lie algebra se(3) to the
Lie group SE(3): given a motion-parameter vector `p = (u; φ) ∈ ℝ⁶`, it
returns the `NodeFrame` `H = (x, φ)` that represents the corresponding rigid
body transformation.

## exp\_SE3 background

The motion-parameter vector `p::VEC6` is partitioned as:

| Components  | Symbol | Meaning                              |
|-------------|--------|--------------------------------------|
| `p[1:3]`    | `u`    | Translational parameter (reduced)    |
| `p[4:6]`    | `φ`    | Cartesian Rotation Vector (CRV)      |

The exponential map formula is:

```
H.phi = φ             (CRV encoding orientation, unchanged)
H.x   = T(-φ) · u    (position vector)
```

where `T(-φ)` is the SO(3) tangent operator ([`T_SO3`](@ref)) evaluated at
`-φ`:

```
T(φ) = I + α(‖φ‖) φ̃/‖φ‖ + β(‖φ‖) (φ̂⊗φ̂ − I)
```

with `α(x) = (cos x − 1)/x` and `β(x) = 1 − sin(x)/x`.

!!! note "Physical interpretation"
    The translational parameter `u` is **not** the position vector directly.
    It is the "reduced translation" that accounts for the coupling between
    rotation and translation in SE(3).  For zero rotation (`φ = 0`), `u`
    equals the position vector exactly.

## exp\_SE3 small-angle limit

Small-angle handling is delegated transparently to [`T_SO3`](@ref):

```
‖φ‖ < 1e-4:  T(-φ) ≈ I − ½φ̃,  so  x ≈ u − ½ φ × u
‖φ‖ = 0  :  T(0) = I,           so  x = u  (pure translation)
```

## exp\_SE3 round-trip

`exp_SE3` and [`log_SE3`](@ref) are exact inverses:

```
exp_SE3(log_SE3(H)) = H   for all H ∈ SE(3)
log_SE3(exp_SE3(p)) = p   for all p ∈ ℝ⁶
```

## exp\_SE3 example

```julia
# Pure translation along x
p = VEC6(1.0, 0.0, 0.0,  0.0, 0.0, 0.0)
H = exp_SE3(p)   # NodeFrame: x = [1,0,0], phi = [0,0,0]

# 90° rotation about z with coupled translation
p = VEC6(1.0, 0.0, 0.0,  0.0, 0.0, π/2)
H = exp_SE3(p)   # H.x = T(-π/2·ẑ)·[1,0,0] ≠ [1,0,0]

# Round-trip
p2 = log_SE3(H)  # recovers p exactly
```

**Dependencies:** [`T_SO3`](@ref)

**See also:** [`log_SE3`](@ref), [`T_SO3`](@ref), [`NodeFrame`](@ref)

```@docs
exp_SE3
```
