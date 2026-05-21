```@meta
CurrentModule = ExpMap
```

# bracket — Anticommutator of Skew-Symmetric Matrices

`bracket` computes the **symmetric (anticommutator) product** of the
skew-symmetric matrices associated with two 3D vectors `u` and `v`.

## bracket formula

```
bracket(u, v) = ũṽ + ṽũ  =  2(u·v) I₃  −  v⊗u  −  u⊗v
```

where `ũ = tilde(u)` and `ṽ = tilde(v)` are the 3×3 skew-symmetric matrices
([`tilde`](@ref)), and `v⊗u` is the outer product ([`outerp`](@ref)).

The result is a **symmetric** 3×3 matrix satisfying `bracket(u, v) = bracket(v, u)`.

## bracket properties

| Property   | Formula                                        |
|------------|------------------------------------------------|
| Symmetric  | `bracket(u, v) = bracket(v, u)`                |
| Bilinear   | `bracket(αu, v) = α · bracket(u, v)`           |
| Trace      | `tr(bracket(u, v)) = 4(u·v)`                  |

## bracket role in SO(3) linearisation

`bracket(u, v)` appears in second-order expansions of SO(3) maps.  In
particular, the Gâteaux derivative of the tangent operator [`T_SO3`](@ref)
with respect to the rotation vector `φ` involves terms of the form
`bracket(k, δk)`, where `k = φ/‖φ‖` is the unit rotation axis and `δk` is
its variation.

**Dependencies:** [`tilde`](@ref)

**See also:** [`tilde`](@ref), [`outerp`](@ref), [`T_SO3`](@ref)

```@docs
bracket
```
