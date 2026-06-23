# ==============================================================================
# MAT6.jl — Immutable 6×6 Matrix Type for Julia
# ==============================================================================
#
# Defines an immutable struct `MAT6` representing a 6×6 matrix stored as six
# column vectors of type `VEC6`.  This column-major storage makes matrix–vector
# and matrix–matrix products particularly efficient.
#
# In rigid-body mechanics, a MAT6 typically represents a spatial operator such
# as a stiffness / compliance / inertia matrix acting on wrenches or twists.
#
# CONSTRUCTION
#   MAT6(c1,…,c6 ::VEC6)               → from six column VEC6 vectors
#   MAT6(v::Vector{Float64})            → from a 36-element flat vector (col-major)
#   MAT6(A::Array{Float64,2})           → from a 6×6 Julia matrix
#   MAT6(A,B,C,D ::MAT3)               → from four 3×3 blocks  [A B; C D]
#   MAT6()                             → zero matrix
#   zero(MAT6)                         → zero matrix
#   one(MAT6)                          → identity matrix  I₆
#
# ARITHMETIC OPERATORS
#   m1 + m2                            → matrix addition
#   m1 - m2                            → matrix subtraction
#   m1 * m2                            → matrix–matrix product
#   m  * v   (VEC6)                    → matrix–vector product  → VEC6
#   s  * m  /  m * s   (Float64)       → scalar multiplication (left or right)
#
# LINEAR ALGEBRA
#   transpose(m)                       → transpose  mᵀ  → MAT6
#   inv(m)                             → matrix inverse  m⁻¹  → MAT6
#   m  \ v   (VEC6)                    → solve  m·x = v  → VEC6
#   m  \ n   (MAT6)                    → solve  m·X = n  → MAT6
#   diag(v::VEC6)                      → diagonal matrix from VEC6 → MAT6
#
# ACCESSORS & INDEXING
#   m.v1 … m.v6                        → j-th column as VEC6
#   m.mat                              → convert to Array{Float64,2}  (6×6)
#   m[j]                               → j-th column  → VEC6
#   m[i, j]                            → scalar element (row i, col j) → Float64
#   setindex(m, val, i, j)             → new MAT6 with element (i,j) set to val
#
# DISPLAY & UTILITIES
#   show(io, m)                        → pretty-print rows in scientific notation
#
# NOTES
#   • MAT6 is immutable: its six VEC6 column fields cannot be mutated in place.
#     `setindex` (without `!`) returns a *new* MAT6 with the updated element.
#   • `inv` and `\` delegate to Julia's built-in dense solver via `.mat`;
#     they require a non-singular matrix.
#   • `diag(v::VEC6)` extends `LinearAlgebra.diag`; requires `using LinearAlgebra`.
#   • `show` uses `@printf`; requires `using Printf`.
#   • This type assumes VEC3, MAT3, VEC6 (see their respective .jl files) are loaded.
# ==============================================================================

using Printf
using LinearAlgebra


# ------------------------------------------------------------------------------
# Struct definition
# ------------------------------------------------------------------------------

"""
    MAT6

Immutable 6×6 matrix stored as six column vectors `v1` … `v6` of type `VEC6`.

Column-major storage makes matrix–vector products (`M * v`) allocation-free and
cache-friendly. Element `(i, j)` is accessed as `M[i, j]` or `M.vj[i]`.
"""
struct MAT6
    v1::VEC6
    v2::VEC6
    v3::VEC6
    v4::VEC6
    v5::VEC6
    v6::VEC6

    # Inner constructor — copies columns defensively to guarantee immutability
    function MAT6(c1::VEC6, c2::VEC6, c3::VEC6, c4::VEC6, c5::VEC6, c6::VEC6)
        new(copy(c1), copy(c2), copy(c3), copy(c4), copy(c5), copy(c6))
    end
end

# Tuple of column field symbols — used by indexing and transpose helpers
const _MAT6_COL_FIELDS  = (:v1, :v2, :v3, :v4, :v5, :v6)
const _VEC6_COMP_FIELDS = (:x1, :x2, :x3, :x4, :x5, :x6)


# ------------------------------------------------------------------------------
# Constructors
# ------------------------------------------------------------------------------

"""
    MAT6() → MAT6

Construct the 6×6 zero matrix.
"""
MAT6() = zero(MAT6)

"""
    zero(MAT6) → MAT6

Return the 6×6 zero matrix (integrates with Julia's `zero` interface).
"""
Base.zero(::Type{MAT6}) = MAT6(ntuple(_ -> zero(VEC6), 6)...)

Base.length(m::MAT6) = 36
Base.broadcastable(m::MAT6) = Ref(m)

"""
    one(MAT6) → MAT6

Return the 6×6 identity matrix I₆ (integrates with Julia's `one` interface).
"""
one(::Type{MAT6}) = MAT6(ntuple(j -> VEC6(ntuple(i -> i == j ? 1.0 : 0.0, 6)), 6)...)

"""
    MAT6(v::Vector{Float64}) → MAT6

Construct a `MAT6` from a 36-element flat `Vector{Float64}` stored in
column-major order (elements 1–6 → column 1, 7–12 → column 2, …).
"""
MAT6(v::Vector{Float64}) = MAT6(
    VEC6(v[ 1: 6]...),
    VEC6(v[ 7:12]...),
    VEC6(v[13:18]...),
    VEC6(v[19:24]...),
    VEC6(v[25:30]...),
    VEC6(v[31:36]...),
)

"""
    MAT6(A::Array{Float64,2}) → MAT6

Construct a `MAT6` from a standard 6×6 Julia matrix.
"""
function MAT6(A::Array{Float64,2})
    @assert size(A) == (6, 6) "Expected a 6×6 matrix, got $(size(A))"
    MAT6(ntuple(j -> VEC6(ntuple(i -> A[i, j], 6)), 6)...)
end

"""
    MAT6(A::MAT3, B::MAT3, C::MAT3, D::MAT3) → MAT6

Construct the block matrix  `[A  B ; C  D]`  from four 3×3 blocks:

    ┌         ┐
    │  A  │ B │   rows 1–3
    │─────┼───│
    │  C  │ D │   rows 4–6
    └         ┘
    cols 1–3    cols 4–6
"""
function MAT6(A::MAT3, B::MAT3, C::MAT3, D::MAT3)
    cols = ntuple(6) do j
        if j <= 3
            # Left block: columns 1–3 from A (top) and C (bottom)
            VEC6(A[1,j], A[2,j], A[3,j], C[1,j], C[2,j], C[3,j])
        else
            # Right block: columns 4–6 from B (top) and D (bottom)
            VEC6(B[1,j-3], B[2,j-3], B[3,j-3], D[1,j-3], D[2,j-3], D[3,j-3])
        end
    end
    MAT6(cols...)
end


# ------------------------------------------------------------------------------
# Arithmetic operators
# ------------------------------------------------------------------------------

"""    m1 + m2  →  MAT6   (column-wise addition)       """
Base.:+(m1::MAT6, m2::MAT6) = MAT6(ntuple(j -> m1[j] + m2[j], 6)...)

"""    m1 - m2  →  MAT6   (column-wise subtraction)    """
Base.:-(m1::MAT6, m2::MAT6) = MAT6(ntuple(j -> m1[j] - m2[j], 6)...)

"""    -m       →  MAT6   (unary negation)              """
Base.:-(m::MAT6) = MAT6(ntuple(j -> -m[j], 6)...)

"""    s * m    →  MAT6   (scalar multiplication, scalar on left)  """
Base.:*(s::Float64, m::MAT6) = MAT6(ntuple(j -> s * getfield(m, j), 6)...)

"""    m * s    →  MAT6   (scalar multiplication, matrix on left)  """
Base.:*(m::MAT6, s::Float64) = s * m

"""
    m * v  →  VEC6   (matrix–vector product)

Computes `(M·v)[i] = Σⱼ M[i,j] · v[j]` for i = 1 … 6.
"""
Base.:*(m::MAT6, v::VEC6) = VEC6(ntuple(i -> sum(m[i,j] * v[j] for j in 1:6), 6)...)

"""
    m1 * m2  →  MAT6   (matrix–matrix product)

Each column of the result is `m1` applied to the corresponding column of `m2`.
"""
Base.:*(m1::MAT6, m2::MAT6) = MAT6(ntuple(j -> m1 * getfield(m2, j), 6)...)

# ------------------------------------------------------------------------------
# Linear algebra
# ------------------------------------------------------------------------------

"""
    transpose(m::MAT6) → MAT6

Return the transpose mᵀ: row i of `m` becomes column i of the result.
"""
function Base.transpose(m::MAT6)
    MAT6(ntuple(i ->
        VEC6(ntuple(j ->
            getfield(getfield(m, _MAT6_COL_FIELDS[j]), _VEC6_COMP_FIELDS[i]),
        6)...),
    6)...)
end

"""
    inv(m::MAT6) → MAT6

Return the inverse m⁻¹ computed via Julia's built-in dense solver.
Requires `m` to be non-singular.
"""
Base.inv(m::MAT6) = MAT6(m.mat \ Matrix{Float64}(I, 6, 6))

"""
    m \\ v  →  VEC6

Solve the linear system `m·x = v` for `x` using Julia's built-in dense solver.
Requires `m` to be non-singular.
"""
function Base.:\(m::MAT6, v::VEC6)
    A = [m[i,j] for i in 1:6, j in 1:6]
    b = [v[i]   for i in 1:6]
    VEC6((A \ b)...)
end

"""
    m \\ n  →  MAT6

Solve the linear matrix system `m·X = n` column by column.
Requires `m` to be non-singular.
"""
function Base.:\(m::MAT6, n::MAT6)
    A = [m[i,j] for i in 1:6, j in 1:6]
    B = [n[i,j] for i in 1:6, j in 1:6]
    MAT6(A \ B)
end

"""
    diag(v::VEC6) → MAT6

Construct the 6×6 diagonal matrix whose diagonal entries are the components of `v`.
Extends `LinearAlgebra.diag`; requires `using LinearAlgebra` at the call site.
"""
function LinearAlgebra.diag(v::VEC6)
    z = 0.0
    MAT6(
        VEC6(v.x1,  z,    z,    z,    z,    z  ),
        VEC6(  z,  v.x2,  z,    z,    z,    z  ),
        VEC6(  z,    z,  v.x3,  z,    z,    z  ),
        VEC6(  z,    z,    z,  v.x4,  z,    z  ),
        VEC6(  z,    z,    z,    z,  v.x5,  z  ),
        VEC6(  z,    z,    z,    z,    z,  v.x6),
    )
end


# ------------------------------------------------------------------------------
# Accessors, indexing, display
# ------------------------------------------------------------------------------

"""
    getproperty(m::MAT6, s::Symbol)

Extend property access so that `m.mat` returns the matrix as an
`Array{Float64,2}` — all other field accesses behave normally.
"""
function Base.getproperty(m::MAT6, s::Symbol)
    if s === :mat
        return [m[i,j] for i in 1:6, j in 1:6]
    end
    return getfield(m, s)
end

"""
    m[j::Int]  →  VEC6

Return the `j`-th column of `m` as a `VEC6` (1-based).
"""
Base.getindex(m::MAT6, j::Int) = getfield(m, j)

"""
    m[i::Int, j::Int]  →  Float64

Return the scalar element at row `i`, column `j` (1-based).
"""
Base.getindex(m::MAT6, i::Int, j::Int) = getfield(getfield(m, j), i)

"""
    setindex(m::MAT6, val::Float64, i::Int, j::Int) → MAT6

Return a new `MAT6` identical to `m` except that element `(i, j)` (1-based)
is replaced by `val`. Since `MAT6` is immutable, the original is unchanged.
"""
function Base.setindex(m::MAT6, val::Float64, i::Int, j::Int)
    new_cols = ntuple(6) do k
        col = getfield(m, k)
        k == j ? VEC6(ntuple(n -> n == i ? val : getfield(col, n), 6)...) : col
    end
    MAT6(new_cols...)
end

"""
    show(io::IO, m::MAT6)

Pretty-print the matrix, one row per line, in scientific notation (`%13.5e`).
Requires `using Printf`.
"""
function Base.show(io::IO, m::MAT6)
    println(io, " ")
    for i in 1:6
        for j in 1:6
            @printf(io, "%13.5e", m[i,j])
        end
        println(io, " ")
    end
end
