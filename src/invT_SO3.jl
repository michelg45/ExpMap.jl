"""
    invT_SO3(φ::RV3) -> MAT3
    invT_SO3(φ::RV3, a::VEC3) -> VEC3

Compute the inverse tangent operator `T⁻¹(φ)` of SO(3), or its action on a vector `a`.

# Background

The inverse tangent operator is the left-trivialized inverse of the tangent map
of the exponential on SO(3).  For a rotation vector `φ ∈ ℝ³` with norm `‖φ‖ = x`,
it reads

    T⁻¹(φ) = I + ½ φ̃ + γ(x) (φ̂ ⊗ φ̂ − I)

where

- `φ̂ = φ / x`                             is the unit rotation axis,
- `φ̃`                                      is the skew-symmetric matrix of `φ`,
- `γ(x) = 1 − (x/2) cot(x/2)`             is the scalar coefficient function
  computed by `invT_functions`.

The relation `T⁻¹(φ) T(φ) = I` holds for all `φ`.

# Small-angle limit

For `‖φ‖ < 1e-4` the first-order approximation

    T⁻¹(φ) ≈ I + ½ φ̃

is returned to avoid numerical cancellation.  The vector form correspondingly
gives `T⁻¹(φ) a ≈ a + ½ (φ × a)`.

# Methods

## `invT_SO3(φ::RV3) -> MAT3`

Return the full 3×3 inverse tangent matrix `T⁻¹(φ)`.

### Arguments
- `φ::RV3`: rotation vector (3-component).

### Returns
- `MAT3`: the inverse tangent matrix `T⁻¹(φ)`.

---

## `invT_SO3(φ::RV3, a::VEC3) -> VEC3`

Return the image `T⁻¹(φ) a` without forming the full matrix.

    T⁻¹(φ) a = a + ½ (φ × a) + γ(x) (φ̂ (φ̂ · a) − a)
             = a + ½ (φ × a) + γ(x)/x² (φ (φ · a) − x² a)

### Arguments
- `φ::RV3`  : rotation vector.
- `a::VEC3` : vector to be mapped.

### Returns
- `VEC3`: the product `T⁻¹(φ) a`.
"""
function invT_SO3(φ::RV3)

    I   = one(MAT3)
    x²  = dotp(φ, φ)           # squared norm ‖φ‖²

    # Small-angle approximation (first-order)
    x² < 1.0e-8 && return I + 0.5 * tilde(φ)

    x   = sqrt(x²)              # norm ‖φ‖
    φ̂   = φ / x                 # unit axis

    γ   = invT_functions(x)     # scalar coefficient

    return I + 0.5 * tilde(φ) + γ * (outerp(φ̂, φ̂) - I)
end

function invT_SO3(φ::RV3, a::VEC3)

    x²  = dotp(φ, φ)           # squared norm ‖φ‖²

    # Small-angle approximation (first-order)
    x² < 1.0e-8 && return a + 0.5 * crossp(φ, a)

    x   = sqrt(x²)              # norm ‖φ‖
    φ   = VEC3(φ)               # ensure concrete type

    γ   = invT_functions(x)     # scalar coefficient

    # Written directly in terms of φ (not φ̂) to avoid recomputing x:
    # γ(x) (φ̂ (φ̂·a) − a)  =  (γ(x)/x²) (φ (φ·a) − x² a)
    return a + 0.5 * crossp(φ, a) + (γ / x²) * (φ * dotp(φ, a) - x² * a)
end
