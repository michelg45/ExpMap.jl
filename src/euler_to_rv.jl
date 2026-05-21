# ==============================================================================
# euler_to_rv.jl — Euler Angles (ZXZ) to Cartesian Rotation Vector
# ==============================================================================
#
# Converts a set of classical ZXZ Euler angles (φ, θ, ψ) into the equivalent
# Cartesian Rotation Vector (CRV) by:
#
#   1. Building the composite rotation matrix via three elementary rotations:
#
#        R(φ, θ, ψ) = R_z(ψ) · R_x(θ) · R_z(φ)
#
#      where R_z(α) is a rotation of angle α about the z-axis (axis 3) and
#      R_x(α) is a rotation of angle α about the x-axis (axis 1).
#
#   2. Extracting the corresponding CRV via the logarithmic map `invR_SO3`.
#
# CONVENTION
#   The ZXZ (or 3-1-3) convention is the classical Euler angle parameterisation
#   used in rigid-body dynamics (Goldstein):
#     φ  — precession angle  (first  rotation about z, axis 3)
#     θ  — nutation  angle   (second rotation about x, axis 1), θ ∈ [0, π]
#     ψ  — spin      angle   (third  rotation about z, axis 3)
#
#   The composition is applied right-to-left: φ first, then θ, then ψ.
#
# SINGULARITY
#   At θ = 0 and θ = π the parameterisation is singular (gimbal lock): only
#   φ + ψ (or φ − ψ) is determined. `invR_SO3` still returns a valid CRV,
#   but the decomposition back into Euler angles is non-unique at these points.
#
# EXTERNAL DEPENDENCIES
#   R_SO3(angle::Float64, axis::Int)  elementary rotation about coordinate axis
#   invR_SO3(R::MAT3)                 logarithmic map SO(3) → so(3)
# ==============================================================================


"""
    euler_to_RV(phi::Float64, theta::Float64, psi::Float64) → RV3

Convert the ZXZ Euler angles `(phi, theta, psi)` to the equivalent Cartesian
Rotation Vector (CRV).

The composite rotation matrix is built as:

    R(φ, θ, ψ) = R_z(ψ) · R_x(θ) · R_z(φ)

where `R_z(α)` and `R_x(α)` are elementary rotations of angle `α` about the
z-axis (axis 3) and x-axis (axis 1), respectively.  The CRV is then recovered
via the logarithmic map [`invR_SO3`](@ref).

## Arguments
- `phi`   : precession angle φ (radians) — first rotation about z
- `theta` : nutation angle θ (radians)   — rotation about x; θ ∈ [0, π]
- `psi`   : spin angle ψ (radians)       — last rotation about z

## Returns
`rv::RV3` — Cartesian Rotation Vector equivalent to R(φ, θ, ψ)

## Example
```julia
julia> rv = euler_to_RV(0.0, π/9, 0.0)
RV3([3.490659e-01, 0.000000e+00, 0.000000e+00])
```

!!! warning "Singularity"
    At θ = 0 and θ = π the ZXZ parameterisation is singular (gimbal lock):
    only φ ± ψ is determined, not φ and ψ individually.  A valid CRV is still
    returned, but the round-trip back to Euler angles is non-unique at these
    points.

**See also:** [`R_SO3`](@ref), [`invR_SO3`](@ref)
"""
function euler_to_RV(phi::Float64, theta::Float64, psi::Float64)
    R = R_SO3(psi, 3) * R_SO3(theta, 1) * R_SO3(phi, 3)
    return invR_SO3(R)
end
