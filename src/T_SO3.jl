"""
    T_SO3(φ::RV3) -> MAT3
    T_SO3(φ::RV3, a::VEC3) -> VEC3

Compute the tangent operator `T(φ)` of SO(3), or its action on a vector `a`.

# Background

The tangent operator (also called the left-trivialized tangent map) arises when
linearising the exponential map on SO(3).  For a rotation vector `φ ∈ ℝ³` with
norm `‖φ‖ = x`, it reads

    T(φ) = I + α(x) φ̂̃ + β(x) (φ̂ ⊗ φ̂ − I)

where

- `φ̂ = φ / x`           is the unit rotation axis,
- `φ̂̃`                   is the skew-symmetric matrix associated with `φ̂`,
- `α(x) = (cos x − 1) / x`
- `β(x) = 1 − sin(x) / x`

are the scalar coefficient functions computed by `T_functions`.

# Small-angle limit

For `‖φ‖ < 1e-4` the first-order approximation

    T(φ) ≈ I − ½ φ̃

is returned to avoid numerical cancellation.  The vector form correspondingly
gives `T(φ) a ≈ a − ½ (φ × a)`.

# Methods

## `T_SO3(φ::RV3) -> MAT3`

Return the full 3×3 tangent matrix `T(φ)`.

### Arguments
- `φ::RV3`: rotation vector (3-component).

### Returns
- `MAT3`: the tangent matrix `T(φ)`.

---

## `T_SO3(φ::RV3, a::VEC3) -> VEC3`

Return the image `T(φ) a` without forming the full matrix.

    T(φ) a = a + α(x) (φ̂ × a) + β(x) (φ̂ (φ̂ · a) − a)

### Arguments
- `φ::RV3`  : rotation vector.
- `a::VEC3` : vector to be mapped.

### Returns
- `VEC3`: the product `T(φ) a`.
"""
function T_SO3(φ::RV3)

    I   = one(MAT3)
    x²  = dotp(φ, φ)           # squared norm ‖φ‖²

    # Small-angle approximation (first-order)
    x² < 1.0e-8 && return I - 0.5 * tilde(φ)

    x   = sqrt(x²)              # norm ‖φ‖
    φ̂   = φ / x                 # unit axis

    α, β = T_functions(x)       # scalar coefficients

    return I + α * tilde(φ̂) + β * (outerp(φ̂, φ̂) - I)
end

function T_SO3(φ::RV3, a::VEC3)

    x²  = dotp(φ, φ)           # squared norm ‖φ‖²

    # Small-angle approximation (first-order)
    x² < 1.0e-8 && return a - 0.5 * crossp(φ, a)

    x   = sqrt(x²)                                  # norm ‖φ‖
    φ̂   = VEC3(φ[1]/x, φ[2]/x, φ[3]/x)             # unit axis

    α, β = T_functions(x)                            # scalar coefficients

    return a + α * crossp(φ̂, a) + β * (φ̂ * dotp(φ̂, a) - a)
end
