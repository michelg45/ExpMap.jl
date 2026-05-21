# ==============================================================================
# MAT3.jl — Immutable 3×3 Matrix Type for Julia
# ==============================================================================
#
# Defines an immutable struct `MAT3` representing a 3×3 matrix with Float64
# entries stored in row-major order (a_ij: row i, column j).
#
# This type is designed to interoperate with `VEC3` (see VEC3.jl) and is
# well-suited for rotation matrices, inertia tensors, and other 3D linear maps.
#
# CONSTRUCTION
#   MAT3(a11,a12,a13, a21,a22,a23, a31,a32,a33)  → from nine Float64 scalars
#   MAT3(m::Array{Float64,2})                     → from a 3×3 Julia matrix
#   MAT3(a::VEC3, b::VEC3, c::VEC3)              → from three column vectors
#   MAT3()                                        → zero matrix
#   zero(MAT3)                                    → zero matrix
#   one(MAT3)                                     → identity matrix  I₃
#
# ARITHMETIC OPERATORS
#   m1 + m2                → matrix addition
#   m1 - m2                → matrix subtraction
#   -m                     → negation (unary minus)
#   m1 * m2                → matrix–matrix product
#   m  * v   (VEC3)        → matrix–vector product  → VEC3
#   m  * s  /  s * m       → scalar multiplication (left or right)
#   m  / s                 → scalar division
#
# LINEAR ALGEBRA
#   transpose(m)           → transpose  mᵀ  → MAT3
#   inv(m)                 → matrix inverse  m⁻¹  → MAT3
#   m \ v  (VEC3)          → solve linear system  m·x = v  → VEC3
#   diag(v::VEC3)          → diagonal matrix from VEC3 components → MAT3
#
# ACCESSORS & INDEXING
#   m.aij                  → direct field access  (e.g. m.a12 → row 1, col 2)
#   m.mat                  → convert to Array{Float64,2}  (3×3 Julia matrix)
#   m[i, j]                → scalar element at row i, column j  → Float64
#   m[:, j]                → j-th column as VEC3
#   setindex(m, val, i, j) → return new MAT3 with element (i,j) set to val
#
# DISPLAY & UTILITIES
#   show(io, m)            → pretty-print rows in scientific notation
#   copy(m)                → deep copy of the matrix
#
# NOTES
#   • MAT3 is immutable: fields cannot be mutated in place.
#     `setindex` returns a *new* MAT3 with the updated element.
#   • `inv` and `\` delegate to Julia's built-in dense linear algebra via
#     the `.mat` accessor; they require a non-singular matrix.
#   • `diag(v::VEC3)` extends `LinearAlgebra.diag` and requires
#     `using LinearAlgebra` at the call site.
#   • This type assumes VEC3 (see VEC3.jl) is already loaded.
#   • `show` uses `@printf` and therefore requires `using Printf`.
# ==============================================================================

using Printf
using LinearAlgebra


# ------------------------------------------------------------------------------
# Struct definition
# ------------------------------------------------------------------------------

"""
    MAT3

Immutable 3×3 matrix with `Float64` entries stored in row-major order.
Fields are named `aij` where `i` is the row index and `j` the column index
(e.g. `a12` is row 1, column 2).
"""
struct MAT3
    a11::Float64;  a12::Float64;  a13::Float64
    a21::Float64;  a22::Float64;  a23::Float64
    a31::Float64;  a32::Float64;  a33::Float64
end

# Flat tuple of field symbols in row-major order — used by indexing helpers
const _MAT3_FIELDS = (
    :a11, :a12, :a13,
    :a21, :a22, :a23,
    :a31, :a32, :a33,
)


# ------------------------------------------------------------------------------
# Constructors
# ------------------------------------------------------------------------------

"""
    MAT3() → MAT3

Construct the 3×3 zero matrix.
"""
MAT3() = MAT3(0.0, 0.0, 0.0,
              0.0, 0.0, 0.0,
              0.0, 0.0, 0.0)

"""
    zero(MAT3) → MAT3

Return the 3×3 zero matrix (integrates with Julia's `zero` interface).
"""
Base.zero(::Type{MAT3}) = MAT3(0.0, 0.0, 0.0,
                                0.0, 0.0, 0.0,
                                0.0, 0.0, 0.0)

"""
    one(MAT3) → MAT3

Return the 3×3 identity matrix I₃ (integrates with Julia's `one` interface).
"""
one(::Type{MAT3}) = MAT3(1.0, 0.0, 0.0,
                          0.0, 1.0, 0.0,
                          0.0, 0.0, 1.0)

"""
    MAT3(m::Array{Float64,2}) → MAT3

Construct a `MAT3` from a standard 3×3 Julia matrix.
"""
MAT3(m::Array{Float64,2}) = MAT3(m[1,1], m[1,2], m[1,3],
                                  m[2,1], m[2,2], m[2,3],
                                  m[3,1], m[3,2], m[3,3])

"""
    MAT3(a::VEC3, b::VEC3, c::VEC3) → MAT3

Construct a `MAT3` whose columns are the three `VEC3` vectors `a`, `b`, `c`.

    ┌ a.x  b.x  c.x ┐
    │ a.y  b.y  c.y │
    └ a.z  b.z  c.z ┘
"""
function MAT3(a::VEC3, b::VEC3, c::VEC3)
    MAT3(a.x, b.x, c.x,
         a.y, b.y, c.y,
         a.z, b.z, c.z)
end


# ------------------------------------------------------------------------------
# Arithmetic operators
# ------------------------------------------------------------------------------

"""    m1 + m2  →  MAT3   (element-wise addition)      """
Base.:+(m1::MAT3, m2::MAT3) = MAT3(
    m1.a11 + m2.a11,  m1.a12 + m2.a12,  m1.a13 + m2.a13,
    m1.a21 + m2.a21,  m1.a22 + m2.a22,  m1.a23 + m2.a23,
    m1.a31 + m2.a31,  m1.a32 + m2.a32,  m1.a33 + m2.a33,
)

"""    m1 - m2  →  MAT3   (element-wise subtraction)   """
Base.:-(m1::MAT3, m2::MAT3) = MAT3(
    m1.a11 - m2.a11,  m1.a12 - m2.a12,  m1.a13 - m2.a13,
    m1.a21 - m2.a21,  m1.a22 - m2.a22,  m1.a23 - m2.a23,
    m1.a31 - m2.a31,  m1.a32 - m2.a32,  m1.a33 - m2.a33,
)

"""    -m       →  MAT3   (unary negation)              """
Base.:-(m::MAT3) = MAT3(
    -m.a11, -m.a12, -m.a13,
    -m.a21, -m.a22, -m.a23,
    -m.a31, -m.a32, -m.a33,
)

"""    m1 * m2  →  MAT3   (matrix–matrix product)       """
Base.:*(m1::MAT3, m2::MAT3) = MAT3(
    m1.a11*m2.a11 + m1.a12*m2.a21 + m1.a13*m2.a31,
    m1.a11*m2.a12 + m1.a12*m2.a22 + m1.a13*m2.a32,
    m1.a11*m2.a13 + m1.a12*m2.a23 + m1.a13*m2.a33,
    m1.a21*m2.a11 + m1.a22*m2.a21 + m1.a23*m2.a31,
    m1.a21*m2.a12 + m1.a22*m2.a22 + m1.a23*m2.a32,
    m1.a21*m2.a13 + m1.a22*m2.a23 + m1.a23*m2.a33,
    m1.a31*m2.a11 + m1.a32*m2.a21 + m1.a33*m2.a31,
    m1.a31*m2.a12 + m1.a32*m2.a22 + m1.a33*m2.a32,
    m1.a31*m2.a13 + m1.a32*m2.a23 + m1.a33*m2.a33,
)

"""    m * v    →  VEC3   (matrix–vector product)       """
Base.:*(m::MAT3, v::VEC3) = VEC3(
    m.a11*v.x + m.a12*v.y + m.a13*v.z,
    m.a21*v.x + m.a22*v.y + m.a23*v.z,
    m.a31*v.x + m.a32*v.y + m.a33*v.z,
)

"""    s * m    →  MAT3   (scalar multiplication, scalar on left)  """
Base.:*(s::Float64, m::MAT3) = MAT3(
    s*m.a11, s*m.a12, s*m.a13,
    s*m.a21, s*m.a22, s*m.a23,
    s*m.a31, s*m.a32, s*m.a33,
)

"""    m * s    →  MAT3   (scalar multiplication, matrix on left)  """
Base.:*(m::MAT3, s::Float64) = s * m

"""    m / s    →  MAT3   (scalar division)                        """
Base.:/(m::MAT3, s::Float64) = MAT3(
    m.a11/s, m.a12/s, m.a13/s,
    m.a21/s, m.a22/s, m.a23/s,
    m.a31/s, m.a32/s, m.a33/s,
)


# ------------------------------------------------------------------------------
# Linear algebra
# ------------------------------------------------------------------------------

"""
    transpose(m::MAT3) → MAT3

Return the transpose mᵀ of the matrix.
"""
Base.transpose(m::MAT3) = MAT3(
    m.a11, m.a21, m.a31,
    m.a12, m.a22, m.a32,
    m.a13, m.a23, m.a33,
)

"""
    inv(m::MAT3) → MAT3

Return the inverse m⁻¹ computed via Julia's built-in dense solver.
Requires `m` to be non-singular.
"""
Base.inv(m::MAT3) = MAT3(m.mat \ [1.0 0.0 0.0;
                                    0.0 1.0 0.0;
                                    0.0 0.0 1.0])

"""
    m \\ v  →  VEC3

Solve the linear system `m·x = v` for `x` using Julia's built-in dense solver.
Requires `m` to be non-singular.
"""
Base.:\(m::MAT3, v::VEC3) = VEC3(m.mat \ v.v)

"""
    diag(v::VEC3) → MAT3

Construct the diagonal matrix whose diagonal entries are the components of `v`:

    diag(v) = diag(v.x, v.y, v.z)

Extends `LinearAlgebra.diag`; requires `using LinearAlgebra` at the call site.
"""
LinearAlgebra.diag(v::VEC3) = MAT3(
    v.x, 0.0, 0.0,
    0.0, v.y, 0.0,
    0.0, 0.0, v.z,
)


# ------------------------------------------------------------------------------
# Accessors, indexing, display, and utilities
# ------------------------------------------------------------------------------

"""
    getproperty(m::MAT3, s::Symbol)

Extend property access so that `m.mat` returns the matrix as an
`Array{Float64,2}` — all other field accesses behave normally.
"""
function Base.getproperty(m::MAT3, s::Symbol)
    s === :mat && return [getfield(m, _MAT3_FIELDS[(i-1)*3 + j]) for i in 1:3, j in 1:3]
    return getfield(m, s)
end

"""
    m[i, j]  →  Float64

Return the scalar element at row `i`, column `j` (1-based).
"""
Base.getindex(m::MAT3, i::Int, j::Int) = getfield(m, _MAT3_FIELDS[(i-1)*3 + j])

"""
    m[:, j]  →  VEC3

Return the `j`-th column of `m` as a `VEC3`.
"""
function Base.getindex(m::MAT3, ::Colon, j::Int)
    VEC3(getfield(m, _MAT3_FIELDS[j]),
         getfield(m, _MAT3_FIELDS[j + 3]),
         getfield(m, _MAT3_FIELDS[j + 6]))
end

"""
    setindex(m::MAT3, val::Float64, i::Int, j::Int) → MAT3

Return a new `MAT3` identical to `m` except that element `(i, j)` (1-based)
is replaced by `val`. Since `MAT3` is immutable, the original is unchanged.
"""
function Base.setindex(m::MAT3, val::Float64, i::Int, j::Int)
    k = 3*(i-1) + j
    fields = ntuple(n -> n == k ? val : getfield(m, _MAT3_FIELDS[n]), 9)
    return MAT3(fields...)
end

"""
    show(io::IO, m::MAT3)

Pretty-print the matrix, one row per line, in scientific notation (`%13.5e`).
Requires `using Printf`.
"""
function Base.show(io::IO, m::MAT3)
    println(io, " ")
    for i in 1:3
        for j in 1:3
            @printf(io, "%13.5e", m[i,j])
        end
        println(io, " ")
    end
end

"""
    copy(m::MAT3) → MAT3

Return a deep copy of `m` (useful when storing `MAT3` inside mutable containers).
"""
Base.copy(m::MAT3) = MAT3(
    m.a11, m.a12, m.a13,
    m.a21, m.a22, m.a23,
    m.a31, m.a32, m.a33,
)
