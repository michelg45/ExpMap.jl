# ==============================================================================
# VEC3.jl — Immutable 3D Vector Type for Julia
# ==============================================================================
#
# Defines a lightweight, immutable struct `VEC3` representing a 3D vector with
# Float64 components, along with a comprehensive set of operations.
#
# CONSTRUCTION
#   VEC3(x, y, z)          → from three Float64 scalars
#   VEC3(v::Vector{Float64})→ from a 3-element Julia vector
#   VEC3()                 → zero vector (0.0, 0.0, 0.0)
#
# ARITHMETIC OPERATORS  (all component-wise unless noted)
#   v1 + v2                → addition
#   v1 - v2                → subtraction
#   -v                     → negation (unary minus)
#   v1 * v2                → component-wise multiplication (Hadamard product)
#   v  * s  /  s * v       → scalar multiplication (left or right)
#   v1 / v2                → component-wise division
#   v  / s                 → scalar division
#   s  / v                 → scalar divided by each component
#
# GEOMETRIC OPERATIONS
#   dotp(v1, v2)           → dot product  → Float64
#   crossp(v1, v2)         → cross product → VEC3
#   norm2(v)               → Euclidean norm (L2) → Float64
#
# ACCESSORS & INDEXING
#   v.x, v.y, v.z          → direct field access
#   v.v                    → convert to Vector{Float64}  [x, y, z]
#   v[1], v[2], v[3]       → integer index (1-based: x, y, z)
#
# UTILITIES
#   copy(v)                → deep copy of the vector
#
# NOTES
#   • VEC3 is immutable: fields cannot be modified after construction.
#   • Component-wise multiplication (*) and division (/) are NOT the dot or
#     cross products; use dotp() and crossp() for those.
#   • All operations preserve Float64 precision.
# ==============================================================================


# ------------------------------------------------------------------------------
# Struct definition
# ------------------------------------------------------------------------------

"""
    VEC3

Immutable 3D vector with `Float64` components `x`, `y`, `z`.
"""
struct VEC3
    x::Float64
    y::Float64
    z::Float64
end


# ------------------------------------------------------------------------------
# Constructors
# ------------------------------------------------------------------------------

"""
    VEC3(a::Vector{Float64}) → VEC3

Construct a `VEC3` from a 3-element `Vector{Float64}`.
"""
VEC3(a::Vector{Float64}) = VEC3(a[1], a[2], a[3])

"""
    VEC3() → VEC3

Construct the zero vector `(0.0, 0.0, 0.0)`.
"""
VEC3() = VEC3(0.0, 0.0, 0.0)


# ------------------------------------------------------------------------------
# Arithmetic operators  (component-wise)
# ------------------------------------------------------------------------------

"""    v1 + v2  →  VEC3   (component-wise addition)     """
Base.:+(v1::VEC3, v2::VEC3) = VEC3(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z)

"""    v1 - v2  →  VEC3   (component-wise subtraction)  """
Base.:-(v1::VEC3, v2::VEC3) = VEC3(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z)

"""    -v       →  VEC3   (unary negation)              """
Base.:-(v::VEC3)             = VEC3(-v.x, -v.y, -v.z)

"""    v1 * v2  →  VEC3   (component-wise / Hadamard product) """
Base.:*(v1::VEC3, v2::VEC3) = VEC3(v1.x * v2.x, v1.y * v2.y, v1.z * v2.z)

"""    v * s    →  VEC3   (scalar multiplication, vector on left)  """
Base.:*(v::VEC3,  s::Float64) = VEC3(v.x * s, v.y * s, v.z * s)

"""    s * v    →  VEC3   (scalar multiplication, scalar on left)  """
Base.:*(s::Float64, v::VEC3)  = VEC3(v.x * s, v.y * s, v.z * s)

"""    v1 / v2  →  VEC3   (component-wise division)     """
Base.:/(v1::VEC3, v2::VEC3)  = VEC3(v1.x / v2.x, v1.y / v2.y, v1.z / v2.z)

"""    v / s    →  VEC3   (scalar division)             """
Base.:/(v::VEC3,  s::Float64) = VEC3(v.x / s, v.y / s, v.z / s)

"""    s / v    →  VEC3   (scalar divided by each component)       """
Base.:/(s::Float64, v::VEC3)  = VEC3(s / v.x, s / v.y, s / v.z)


# ------------------------------------------------------------------------------
# Geometric operations
# ------------------------------------------------------------------------------

"""
    dotp(v1::VEC3, v2::VEC3) → Float64

Return the dot product  `v1 · v2 = v1.x·v2.x + v1.y·v2.y + v1.z·v2.z`.
"""
dotp(v1::VEC3, v2::VEC3) = v1.x * v2.x + v1.y * v2.y + v1.z * v2.z

"""
    crossp(v1::VEC3, v2::VEC3) → VEC3

Return the cross product `v1 × v2` as a new `VEC3`.
"""
crossp(v1::VEC3, v2::VEC3) = VEC3(
    v1.y * v2.z - v1.z * v2.y,
    v1.z * v2.x - v1.x * v2.z,
    v1.x * v2.y - v1.y * v2.x,
)

"""
    norm2(v::VEC3) → Float64

Return the Euclidean (L2) norm  `‖v‖ = √(x² + y² + z²)`.
"""
norm2(v::VEC3) = sqrt(v.x^2 + v.y^2 + v.z^2)


# ------------------------------------------------------------------------------
# Accessors, indexing, and utilities
# ------------------------------------------------------------------------------

# Mapping from integer index to field symbol (1-based)
const _VEC3_FIELDS = (:x, :y, :z)

"""
    getproperty(v::VEC3, s::Symbol)

Extend property access so that `v.v` returns the vector as a
`Vector{Float64}` — all other field accesses behave normally.
"""
function Base.getproperty(v::VEC3, s::Symbol)
    s === :v && return [v.x, v.y, v.z]
    return getfield(v, s)
end

"""
    v[i::Int]  →  Float64

1-based integer indexing: `v[1]` → x, `v[2]` → y, `v[3]` → z.
"""
Base.getindex(v::VEC3, i::Int) = getfield(v, _VEC3_FIELDS[i])

"""
    copy(v::VEC3) → VEC3

Return a deep copy of `v` (useful when storing `VEC3` inside mutable containers).
"""
Base.copy(v::VEC3) = VEC3(copy(v.x), copy(v.y), copy(v.z))
