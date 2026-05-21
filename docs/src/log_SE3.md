```@meta
CurrentModule = ExpMap
```

# log_SE3 — Logarithmic Map SE(3)

`log_SE3` implements the **logarithmic map** from the Lie group SE(3) back to
its Lie algebra se(3): given a `NodeFrame` `H = (x, φ)`, it recovers the
motion-parameter vector `p = (u; φ) ∈ ℝ⁶`.

## log\_SE3 background

The logarithmic map formula is:

```
p[4:6] = H.phi                (CRV, unchanged)
p[1:3] = T⁻¹(-H.phi) · H.x  (translational parameter)
```

where `T⁻¹(-φ)` is the inverse of the SO(3) tangent operator
([`invT_SO3`](@ref)) evaluated at `-φ`:

```
T⁻¹(φ) = I + ½ φ̃ + γ(‖φ‖) (φ̂⊗φ̂ − I)
```

with `γ(x) = 1 − (x/2) cot(x/2)`.

The inverse tangent operator `T⁻¹(-φ)` recovers the translational parameter
`u` from the position vector `x`:

```
u = T⁻¹(-φ) · x   ⟺   x = T(-φ) · u
```

## log\_SE3 small-angle limit

Small-angle handling is delegated transparently to [`invT_SO3`](@ref):

```
‖φ‖ < 1e-4:  T⁻¹(-φ) ≈ I + ½φ̃,  so  u ≈ x + ½ φ × x
‖φ‖ = 0  :  T⁻¹(0) = I,           so  u = x  (pure translation)
```

## log\_SE3 round-trip

`log_SE3` and [`exp_SE3`](@ref) are exact inverses:

```
log_SE3(exp_SE3(p)) = p   for all p ∈ ℝ⁶
exp_SE3(log_SE3(H)) = H   for all H ∈ SE(3)
```

## log\_SE3 example

```julia
# Build a frame and recover its motion parameters
x   = VEC3(1.0, 2.0, 0.0)
phi = RV3(0.0, 0.0, π/4)
H   = NodeFrame(x, phi)

p = log_SE3(H)    # VEC6: p[4:6] = phi,  p[1:3] = T⁻¹(-phi)·x

# Round-trip
H2 = exp_SE3(p)   # recovers H exactly
```

**Dependencies:** [`invT_SO3`](@ref)

**See also:** [`exp_SE3`](@ref), [`invT_SO3`](@ref), [`NodeFrame`](@ref)

```@docs
log_SE3
```
