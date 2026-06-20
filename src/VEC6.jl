# ==============================================================================
# VEC6.jl — Immutable 6-Component Vector Type for Julia
# ==============================================================================
#
# Defines an immutable struct `VEC6` representing a 6-component vector with
# Float64 entries, indexed x1 … x6.
#
# In the context of rigid-body mechanics, a VEC6 typically encodes a spatial
# (6D) quantity such as:
#   • a wrench  : [force(3); torque(3)]
#   • a twist   : [linear velocity(3); angular velocity(3)]
#   • a state   : [position(3); rotation vector(3)]
#
# CONSTRUCTION
#   VEC6(x1,x2,x3,x4,x5,x6)    → from six Float64 scalars
#   VEC6(t::NTuple{6,Float64})  → from a 6-tuple
#   VEC6(v::Vector{Float64})    → from a 6-element Julia vector (asserts length)
#   VEC6(a::VEC3, b::VEC3)      → concatenate two VEC3   [a; b]
#   VEC6(a::VEC3, b::RV3)       → concatenate VEC3 and RV3  [a; b]
#   zero(VEC6)                  → zero vector  (0,0,0,0,0,0)
#
# ARITHMETIC OPERATORS  (all component-wise)
#   a + b                       → addition
#   a - b                       → subtraction
#   -v                          → negation (unary minus)
#   λ * v  /  v * λ             → scalar multiplication (Real, left or right)
#   v / λ                       → scalar division
#
# GEOMETRIC OPERATIONS
#   norm2(v)                    → Euclidean norm  ‖v‖  → Float64
#
# UTILITY CONSTRUCTORS
#   VEC6_unit(k, a)             → vector with component k set to a, rest zero
#
# ACCESSORS & INDEXING
#   v.x1 … v.x6                → direct field access
#   v.v                         → convert to Vector{Float64}  [x1,…,x6]
#   v[i]                        → 1-based integer indexing  → Float64
#   setindex(v, val, i)         → return new VEC6 with component i set to val
#
# DISPLAY & UTILITIES
#   show(io, v)                 → prints  VEC6(x1, x2, x3, x4, x5, x6)
#   copy(v)                     → deep copy of the vector
#
# NOTES
#   • VEC6 is immutable: fields cannot be mutated in place.
#     `setindex` (without `!`) returns a *new* VEC6 with the updated component.
#   • Scalar multiplication accepts any `Real` (not just `Float64`), which
#     allows convenient use of integer literals (e.g. `2 * v`).
#   • This type assumes VEC3 and RV3 (see VEC3.jl, RV3.jl) are already loaded.
# ==============================================================================


# ------------------------------------------------------------------------------
# Struct definition
# ------------------------------------------------------------------------------

"""
    VEC6

Immutable 6-component vector with `Float64` entries `x1` … `x6`.

Commonly used to represent spatial (6D) quantities in rigid-body mechanics,
such as wrenches [force; torque] or twists [linear vel.; angular vel.].
"""
struct VEC6
    x1::Float64
    x2::Float64
    x3::Float64
    x4::Float64
    x5::Float64
    x6::Float64
end


# ------------------------------------------------------------------------------
# Constructors
# ------------------------------------------------------------------------------

"""
    VEC6(t::NTuple{6,Float64}) → VEC6

Construct a `VEC6` from a 6-tuple of `Float64`.
"""
VEC6(t::NTuple{6, Float64}) = VEC6(t...)

"""
    VEC6(v::Vector{Float64}) → VEC6

Construct a `VEC6` from a 6-element `Vector{Float64}`.
Asserts that the vector has exactly 6 elements.
"""
VEC6(v::Vector{Float64}) = (@assert length(v) == 6; VEC6(v...))

"""
    VEC6(a::VEC3, b::VEC3) → VEC6

Concatenate two `VEC3` into a `VEC6`: `[a.x, a.y, a.z, b.x, b.y, b.z]`.
"""
VEC6(a::VEC3, b::VEC3) = VEC6(a.x, a.y, a.z, b.x, b.y, b.z)

"""
    VEC6(a::VEC3, b::RV3) → VEC6

Concatenate a `VEC3` and an `RV3` into a `VEC6`: `[a.x, a.y, a.z, b.x, b.y, b.z]`.
"""
VEC6(a::VEC3, b::RV3) = VEC6(a.x, a.y, a.z, b.x, b.y, b.z)

"""
    zero(VEC6) → VEC6

Return the 6-component zero vector (integrates with Julia's `zero` interface).
"""
Base.zero(::Type{VEC6}) = VEC6(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

Base.length(a::VEC6) = 6
Base.broadcastable(v::VEC6) = Ref(v)

# ------------------------------------------------------------------------------
# Arithmetic operators  (component-wise)
# ------------------------------------------------------------------------------

"""    a + b  →  VEC6   (component-wise addition)       """
Base.:+(a::VEC6, b::VEC6) = VEC6(ntuple(i ->  a[i] + b[i], 6)...)

"""    a - b  →  VEC6   (component-wise subtraction)    """
Base.:-(a::VEC6, b::VEC6) = VEC6(ntuple(i ->  a[i] - b[i], 6)...)

"""    -v     →  VEC6   (unary negation)                """
Base.:-(v::VEC6)           = VEC6(ntuple(i -> -v[i],        6)...)

"""    λ * v  →  VEC6   (scalar multiplication, scalar on left;  λ::Real)  """
Base.:*(λ::Real, v::VEC6)  = VEC6(ntuple(i ->  λ * v[i],   6)...)

"""    v * λ  →  VEC6   (scalar multiplication, vector on left;  λ::Real)  """
Base.:*(v::VEC6, λ::Real)  = λ * v

"""    v / λ  →  VEC6   (scalar division;  λ::Real)     """
Base.:/(v::VEC6, λ::Real)  = VEC6(ntuple(i ->  v[i] / λ,   6)...)


# ------------------------------------------------------------------------------
# Geometric operations
# ------------------------------------------------------------------------------

"""
    norm2(v::VEC6) → Float64

Return the Euclidean norm  `‖v‖ = √(x1² + x2² + x3² + x4² + x5² + x6²)`.
"""
norm2(v::VEC6) = sqrt(v.x1^2 + v.x2^2 + v.x3^2 + v.x4^2 + v.x5^2 + v.x6^2)


# ------------------------------------------------------------------------------
# Utility constructors
# ------------------------------------------------------------------------------

"""
    VEC6_unit(k::Int, a::Float64) → VEC6

Return a `VEC6` with component `k` set to `a` and all other components zero.
Useful for constructing canonical basis vectors scaled by `a`.
"""
function VEC6_unit(k::Int, a::Float64)
    VEC6(ntuple(i -> i == k ? a : 0.0, 6)...)
end


# ------------------------------------------------------------------------------
# Accessors, indexing, display, and utilities
# ------------------------------------------------------------------------------

"""
    getproperty(v::VEC6, s::Symbol)

Extend property access so that `v.v` returns the vector as a
`Vector{Float64}` — all other field accesses behave normally.
"""
function Base.getproperty(v::VEC6, s::Symbol)
    s === :v && return [getfield(v, i) for i in 1:6]
    return getfield(v, s)
end

"""
    v[i::Int]  →  Float64

1-based integer indexing: `v[1]` → x1, …, `v[6]` → x6.
"""
Base.getindex(v::VEC6, i::Int) = getfield(v, i)

"""
    setindex(v::VEC6, val::Float64, i::Int) → VEC6

Return a new `VEC6` identical to `v` except that component `i` (1-based)
is replaced by `val`. Since `VEC6` is immutable, the original is unchanged.
"""
function Base.setindex(v::VEC6, val::Float64, i::Int)
    fields = ntuple(k -> k == i ? val : getfield(v, k), 6)
    return VEC6(fields...)
end

"""
    show(io::IO, v::VEC6)

Pretty-print a `VEC6` as `VEC6(x1, x2, x3, x4, x5, x6)`.
"""
Base.show(io::IO, v::VEC6) = print(io, "VEC6(", join([v[i] for i in 1:6], ", "), ")")

"""
    copy(v::VEC6) → VEC6

Return a deep copy of `v` (useful when storing `VEC6` inside mutable containers).
"""
Base.copy(v::VEC6) = VEC6(v.x1, v.x2, v.x3, v.x4, v.x5, v.x6)
