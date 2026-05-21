# ==============================================================================
# bracket.jl — Anticommutator of Skew-Symmetric Matrices (SO(3) bracket)
# ==============================================================================
#
# Computes the symmetric (anticommutator) product of two skew-symmetric 3×3
# matrices associated with 3D vectors u and v:
#
#   bracket(u, v) = ũṽ + ṽũ
#
# Using the identity  (ũṽ)ᵢⱼ = (u·v) δᵢⱼ − vᵢ uⱼ,  this reduces to:
#
#   bracket(u, v) = 2(u·v) I₃  −  v⊗u  −  u⊗v
#
# The result is a symmetric 3×3 matrix (since bracket(u,v) = bracket(v,u)).
#
# PROPERTIES
#   • Symmetric:   bracket(u, v) = bracket(v, u)ᵀ = bracket(v, u)
#   • Bilinear:    bracket(αu, v) = α · bracket(u, v)
#   • trace:       tr(bracket(u, v)) = 2(u·v) tr(I) − 2(u·v) = 4(u·v)
#
# ROLE IN SO(3) LINEARISATION
#
#   bracket(u, v) appears in the second-order expansion of SO(3) maps,
#   for example in the linearisation of the tangent operator T_SO3 with
#   respect to the rotation vector φ:
#
#   The Gâteaux derivative of T(φ) in direction δφ involves terms of the
#   form bracket(k, δk) where k = φ/‖φ‖ is the unit rotation axis.
#
# METHOD
#   bracket(u::VEC3, v::VEC3) → MAT3
#
# EXTERNAL DEPENDENCIES
#   tilde(v::VEC3) → MAT3   skew-symmetric matrix
#   VEC3, MAT3               see VEC3.jl, MAT3.jl
#
# ==============================================================================


"""
    bracket(u::VEC3, v::VEC3) → MAT3

Return the **anticommutator** (symmetric product) of the skew-symmetric
matrices of `u` and `v`:

    bracket(u, v) = ũṽ + ṽũ  =  2(u·v) I₃ − v⊗u − u⊗v

The result is a **symmetric** 3×3 matrix.  The identity
`bracket(u, v) = bracket(v, u)` holds for all `u, v ∈ ℝ³`.

This product appears in second-order expansions of SO(3) maps,
notably in the Gâteaux derivative of the tangent operator [`T_SO3`](@ref).
"""
bracket(u::VEC3, v::VEC3) = tilde(u) * tilde(v) + tilde(v) * tilde(u)
