# ==============================================================================
# mat6_array.jl — Conversion of MAT6 Collections to Plain Float64 Arrays
# ==============================================================================
#
# Provides two conversion functions that unpack collections of `MAT6` objects
# into standard Julia `Array{Float64,2}` (matrix) by assembling the individual
# 6×6 blocks into a single dense block matrix.
#
# LAYOUT CONVENTIONS
#   vec_mat6  :  Vector{MAT6}  (length N)   →  Array{Float64,2}  (6 × 6N)
#                V[k] occupies columns  (k-1)·6+1 : k·6  of the result.
#
#   mat_mat6  :  Matrix{MAT6}  (size M×N)   →  Array{Float64,2}  (6M × 6N)
#                M[p,q] occupies rows  (p-1)·6+1 : p·6
#                         and columns  (q-1)·6+1 : q·6  of the result.
#
# In both cases the 6×6 block layout is preserved, so the result can be used
# directly as a conventional dense matrix (I/O, linear solvers, plotting).
#
# EXTERNAL DEPENDENCIES
#   MAT6 (see MAT6.jl)  —  provides the getindex(m, i, j) accessor
# ==============================================================================


"""
    vec_mat6(V::Vector{MAT6}) → Array{Float64,2}

Assemble a length-`N` vector of `MAT6` matrices into a single
`6 × 6N` block-row matrix.

Block `k` (1-based) of the result contains `V[k]`:
- rows    : 1:6
- columns : (k-1)·6+1 : k·6

## Arguments
- `V` : vector of `N` `MAT6` matrices

## Returns
`A::Array{Float64,2}` of size `(6, 6N)`

## Example
```julia
julia> V = [one(MAT6), zero(MAT6)];
julia> A = vec_mat6(V);          # size (6, 12)
julia> A[:, 1:6] == one(MAT6).mat
true
```

**See also:** [`mat_mat6`](@ref), [`MAT6`](@ref)
"""
function vec_mat6(V::Vector{MAT6})
    N = length(V)
    A = Array{Float64,2}(undef, 6, 6*N)
    for k in 1:N
        c0 = (k-1)*6
        for j in 1:6
            for i in 1:6
                A[i, c0+j] = V[k][i, j]
            end
        end
    end
    return A
end


"""
    mat_mat6(M::Matrix{MAT6}) → Array{Float64,2}

Assemble an `(m × n)` matrix of `MAT6` blocks into a single
`6m × 6n` dense block matrix.

Block `(p, q)` (1-based) of the result contains `M[p, q]`:
- rows    : (p-1)·6+1 : p·6
- columns : (q-1)·6+1 : q·6

## Arguments
- `M` : `m × n` matrix of `MAT6` objects

## Returns
`A::Array{Float64,2}` of size `(6m, 6n)`

## Example
```julia
julia> M = [one(MAT6) zero(MAT6); zero(MAT6) one(MAT6)];   # 2×2 block-identity
julia> A = mat_mat6(M);          # size (12, 12)
julia> A[1:6, 1:6] == one(MAT6).mat
true
```

**See also:** [`vec_mat6`](@ref), [`MAT6`](@ref)
"""
function mat_mat6(M::Matrix{MAT6})
    (m, n) = size(M)
    A = Array{Float64,2}(undef, 6*m, 6*n)
    for q in 1:n
        c0 = (q-1)*6
        for p in 1:m
            r0 = (p-1)*6
            for j in 1:6
                for i in 1:6
                    A[r0+i, c0+j] = M[p, q][i, j]
                end
            end
        end
    end
    return A
end
