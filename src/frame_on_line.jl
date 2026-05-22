# ==============================================================================
# frame_on_line.jl — Orthonormal Frame Aligned on a Line Segment
# ==============================================================================
#
# Constructs an orthonormal frame (returned as a CRV) whose first axis n_1
# is aligned along the line segment (x_1, x_2).  The two remaining axes
# n_2 and n_3 are built by a randomised Gram-Schmidt procedure:
#
#   1. n_1 = (x_2 − x_1) / ‖x_2 − x_1‖
#   2. Choose a random unit vector r not nearly parallel to n_1
#      (|n_1 · r| ≤ 0.5, up to 5 attempts).
#   3. n_3 = normalise(n_1 × r)
#   4. n_2 = normalise(n_3 × n_1)
#   5. Return invR_SO3([n_1 | n_2 | n_3]).
#
# NOTES
#   • n_2 and n_3 are not deterministic: they depend on the random draw.
#     If a reproducible frame is required, seed the RNG before calling
#     (e.g. Random.seed!(n)).
#   • The algorithm is robust for any non-degenerate segment (x_1 ≠ x_2).
#
# EXTERNAL DEPENDENCIES
#   VEC3, norm2, dotp, crossp   — see VEC3.jl
#   MAT3                        — see MAT3.jl
#   invR_SO3                    — logarithmic map SO(3) → so(3)
# ==============================================================================


"""
    frame_on_line(x_1::VEC3, x_2::VEC3) → RV3

Build an orthonormal frame whose first axis **n₁** is aligned along the line
segment from `x_1` to `x_2`, and return its Cartesian Rotation Vector.

The two complementary axes **n₂**, **n₃** are determined by a randomised
Gram-Schmidt procedure:

1. `n₁ = (x₂ − x₁) / ‖x₂ − x₁‖`
2. Pick a random unit vector `r` with `|n₁ · r| ≤ 0.5` (up to 5 attempts).
3. `n₃ = normalise(n₁ × r)`
4. `n₂ = normalise(n₃ × n₁)`
5. Return [`invR_SO3`](@ref)`([n₁ | n₂ | n₃])`.

## Arguments
- `x_1`, `x_2` : endpoint positions of the line segment (`VEC3`)

## Returns
`rv::RV3` — rotation vector of the orthonormal frame aligned on the line

!!! note "Non-deterministic output"
    **n₂** and **n₃** depend on the random draw.  For a reproducible frame,
    seed Julia's RNG with `Random.seed!(n)` before calling.

**See also:** [`invR_SO3`](@ref), [`R_SO3`](@ref), [`NodeFrame`](@ref)
"""
function frame_on_line(x_1::VEC3, x_2::VEC3)
    dx  = x_2 - x_1
    n_1 = dx / norm2(dx)

    n_2 = VEC3()
    for _ in 1:5
        rd  = VEC3(rand(3))
        n_2 = rd / norm2(rd)
        abs(dotp(n_1, n_2)) <= 0.5 && break
    end

    n_3 = crossp(n_1, n_2)
    n_3 = n_3 / norm2(n_3)
    n_2 = crossp(n_3, n_1)
    n_2 = n_2 / norm2(n_2)

    return invR_SO3(MAT3(n_1, n_2, n_3))
end
