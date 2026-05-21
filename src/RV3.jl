# ==============================================================================
# RV3.jl — Cartesian Rotation Vector (CRV) Type for Julia
# ==============================================================================
#
# Defines an immutable struct `RV3` representing a Cartesian Rotation Vector
# (CRV), also known as a rotation pseudo-vector or Euler vector.
#
# A CRV encodes a 3D rotation as  ψ = θ · n̂,  where:
#   • θ  is the rotation angle (‖ψ‖, in radians)
#   • n̂  is the unit rotation axis  (ψ / ‖ψ‖)
#
# The zero vector RV3() represents the identity (no rotation).
#
# CONSTRUCTION
#   RV3(x, y, z)             → from three Float64 scalars
#   RV3(v::Vector{Float64})  → from a 3-element Julia vector
#   RV3(v::VEC3)             → from a VEC3 (reinterpret as CRV)
#   RV3()                    → identity rotation  (0.0, 0.0, 0.0)
#
# CROSS-TYPE CONVERSION
#   VEC3(r::RV3)             → reinterpret CRV as a plain VEC3
#   RV3(v::VEC3)             → reinterpret VEC3 as a CRV
#
# SCALAR OPERATORS
#   r * s  /  s * r          → scalar multiplication (left or right)
#   r / s                    → scalar division
#   -r                       → negation (reverses rotation direction)
#
# GEOMETRIC OPERATIONS
#   dotp(r1, r2)             → dot product  (RV3 · RV3)  → Float64
#   dotp(r,  v)              → dot product  (RV3 · VEC3) → Float64
#   crossp(r1, r2)           → cross product (RV3 × RV3) → RV3
#   crossp(r,  v)            → cross product (RV3 × VEC3)→ VEC3
#   norm2(r)                 → rotation angle ‖r‖ (Euclidean norm) → Float64
#
# ACCESSORS & INDEXING
#   r.x, r.y, r.z            → direct field access
#   r.v                      → convert to Vector{Float64}  [x, y, z]
#   r[1], r[2], r[3]         → 1-based integer indexing (x, y, z)
#   setindex(r, val, i)      → return new RV3 with component i set to val
#
# DISPLAY
#   show(io, r)              → prints  RV3(x, y, z)
#
# UTILITIES
#   copy(r)                  → deep copy of the rotation vector
#
# NOTES
#   • RV3 is immutable: components cannot be mutated in place.
#     `setindex` returns a *new* RV3 with the updated component.
#   • Scalar addition/subtraction and component-wise * or / between two RV3
#     values are intentionally not defined: such operations are geometrically
#     meaningless for rotation vectors.
#   • This type assumes VEC3 (see VEC3.jl) is already loaded.
# ==============================================================================


# ------------------------------------------------------------------------------
# Struct definition
# ------------------------------------------------------------------------------

"""
    RV3

Immutable Cartesian Rotation Vector (CRV) with `Float64` components `x`, `y`, `z`.

The vector encodes a 3D rotation as  `ψ = θ · n̂`, where `‖ψ‖` is the rotation
angle in radians and `ψ / ‖ψ‖` is the unit rotation axis.
The zero vector `RV3()` represents the identity (no rotation).
"""
struct RV3
    x::Float64
    y::Float64
    z::Float64
end


# ------------------------------------------------------------------------------
# Constructors
# ------------------------------------------------------------------------------

"""
    RV3() → RV3

Construct the identity rotation vector `(0.0, 0.0, 0.0)`.
"""
RV3() = RV3(0.0, 0.0, 0.0)

"""
    RV3(a::Vector{Float64}) → RV3

Construct an `RV3` from a 3-element `Vector{Float64}`.
"""
RV3(a::Vector{Float64}) = RV3(a[1], a[2], a[3])

"""
    RV3(a::VEC3) → RV3

Reinterpret a `VEC3` as an `RV3` (component-wise copy, no geometric transform).
"""
RV3(a::VEC3) = RV3(a.x, a.y, a.z)

"""
    VEC3(a::RV3) → VEC3

Reinterpret an `RV3` as a `VEC3` (component-wise copy, no geometric transform).
"""
VEC3(a::RV3) = VEC3(a.x, a.y, a.z)


# ------------------------------------------------------------------------------
# Scalar operators
# ------------------------------------------------------------------------------

"""    r * s  →  RV3   (scalar multiplication, vector on left)   """
Base.:*(v::RV3, s::Float64)  = RV3(v.x * s, v.y * s, v.z * s)

"""    s * r  →  RV3   (scalar multiplication, scalar on left)   """
Base.:*(s::Float64, v::RV3)  = RV3(v.x * s, v.y * s, v.z * s)

"""    r / s  →  RV3   (scalar division)                         """
Base.:/(v::RV3, s::Float64)  = RV3(v.x / s, v.y / s, v.z / s)

"""    -r     →  RV3   (negation — reverses the rotation direction) """
Base.:-(v::RV3)              = RV3(-v.x, -v.y, -v.z)


# ------------------------------------------------------------------------------
# Geometric operations
# ------------------------------------------------------------------------------

"""
    dotp(r1::RV3, r2::RV3) → Float64

Dot product of two rotation vectors: `r1 · r2`.
"""
dotp(v1::RV3, v2::RV3)  = v1.x * v2.x + v1.y * v2.y + v1.z * v2.z

"""
    dotp(r::RV3, v::VEC3) → Float64

Dot product of a rotation vector with a plain vector: `r · v`.
"""
dotp(v1::RV3, v2::VEC3) = v1.x * v2.x + v1.y * v2.y + v1.z * v2.z

"""
    crossp(r1::RV3, r2::RV3) → RV3

Cross product of two rotation vectors: `r1 × r2`.
"""
crossp(v1::RV3, v2::RV3) = RV3(
    v1.y * v2.z - v1.z * v2.y,
    v1.z * v2.x - v1.x * v2.z,
    v1.x * v2.y - v1.y * v2.x,
)

"""
    crossp(r::RV3, v::VEC3) → VEC3

Cross product of a rotation vector with a plain vector: `r × v`.
Returns a `VEC3` since the result is a plain spatial vector.
"""
crossp(v1::RV3, v2::VEC3) = VEC3(
    v1.y * v2.z - v1.z * v2.y,
    v1.z * v2.x - v1.x * v2.z,
    v1.x * v2.y - v1.y * v2.x,
)

"""
    norm2(r::RV3) → Float64

Return the Euclidean norm `‖r‖ = √(x² + y² + z²)`, which equals the
rotation angle θ (in radians).
"""
norm2(v::RV3) = sqrt(v.x^2 + v.y^2 + v.z^2)


# ------------------------------------------------------------------------------
# Accessors, indexing, display, and utilities
# ------------------------------------------------------------------------------

# Mapping from 1-based integer index to field symbol
const _RV3_FIELDS = (:x, :y, :z)

"""
    getproperty(r::RV3, s::Symbol)

Extend property access so that `r.v` returns the vector as a
`Vector{Float64}` — all other field accesses behave normally.
"""
function Base.getproperty(a::RV3, s::Symbol)
    s === :v && return [a.x, a.y, a.z]
    return getfield(a, s)
end

"""
    r[i::Int]  →  Float64

1-based integer indexing: `r[1]` → x, `r[2]` → y, `r[3]` → z.
"""
Base.getindex(v::RV3, i::Int) = getfield(v, _RV3_FIELDS[i])

"""
    setindex(r::RV3, val::Float64, i::Int) → RV3

Return a new `RV3` identical to `r` except that component `i` (1-based)
is replaced by `val`. Since `RV3` is immutable, the original is unchanged.
"""
function Base.setindex(a::RV3, v::Float64, i::Int)
    fields = ntuple(k -> k == i ? v : getfield(a, _RV3_FIELDS[k]), 3)
    return RV3(fields...)
end

"""
    show(io::IO, r::RV3)

Pretty-print an `RV3` as `RV3(x, y, z)`.
"""
function Base.show(io::IO, v::RV3)
    print(io, "RV3(", v.x, ", ", v.y, ", ", v.z, ")")
end

"""
    copy(r::RV3) → RV3

Return a deep copy of `r` (useful when storing `RV3` inside mutable containers).
"""
Base.copy(r::RV3) = RV3(copy(r.x), copy(r.y), copy(r.z))
