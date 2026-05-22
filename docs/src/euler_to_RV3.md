```@meta
CurrentModule = ExpMap
```

# euler\_to\_RV3 — Euler Angles (ZXZ) to CRV

`euler_to_RV3` converts a set of **classical ZXZ Euler angles** (φ, θ, ψ) to
the equivalent **Cartesian Rotation Vector** (CRV) ψ ∈ ℝ³.

## ZXZ Euler angle convention

The ZXZ (or 3-1-3) convention is the classical Euler angle parameterisation
used in rigid-body dynamics (Goldstein).  Three successive elementary rotations
are composed:

```
R(φ, θ, ψ) = R_z(ψ) · R_x(θ) · R_z(φ)
```

| Angle | Name       | Axis  | Range    |
|-------|------------|-------|----------|
| φ     | precession | z (3) | [0, 2π)  |
| θ     | nutation   | x (1) | [0, π]   |
| ψ     | spin       | z (3) | [0, 2π)  |

The composition is applied right-to-left: φ first (about z), then θ (about x),
then ψ (about z again).

## Algorithm

1. Build the composite rotation matrix R(φ, θ, ψ) = R_z(ψ) · R_x(θ) · R_z(φ)
   using [`R_SO3`](@ref) for each elementary rotation.
2. Extract the CRV via the logarithmic map [`invR_SO3`](@ref).

## Singularity (gimbal lock)

At **θ = 0** and **θ = π** the ZXZ parameterisation is singular: only φ + ψ
(or φ − ψ) is determined, not φ and ψ individually.  `euler_to_RV3` still
returns a valid CRV via `invR_SO3`, but the decomposition back to Euler angles
is non-unique at these points.

**Dependencies:** [`R_SO3`](@ref), [`invR_SO3`](@ref)

**See also:** [`R_SO3`](@ref), [`invR_SO3`](@ref)

```@docs
euler_to_RV3
```
