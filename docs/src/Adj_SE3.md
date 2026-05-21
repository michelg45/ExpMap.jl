```@meta
CurrentModule = ExpMap
```

# Adj_SE3 — Adjoint Representation of SE(3)

`Adj_SE3` provides four overloads implementing the group adjoint `Ad(H)`,
the matrix-free adjoint action `Ad(H)·a`, the small adjoint `ad(a)`, and
the Lie bracket `[a, b]` of se(3).

## Adj\_SE3 group adjoint

For `H = (x, φ) ∈ SE(3)` with rotation matrix `R = R_SO3(φ)`:

```
Ad(H) = [R    x̃·R]   ∈ ℝ⁶ˣ⁶
        [0    R  ]
```

where `x̃ = tilde(x)`.  `Ad(H)` changes the reference frame of a twist: if
`a` is expressed in body frame, `Ad(H)·a` is the same twist in spatial frame.

The matrix-free form (second overload) avoids allocating the full 6×6 matrix:

```
(Ad(H)·a)[1:3] = R·a_u  +  x × (R·a_ω)
(Ad(H)·a)[4:6] = R·a_ω
```

!!! tip "Efficiency"
    Prefer `Adj_SE3(H, a)` over `Adj_SE3(H) * a` when only the image vector
    is needed.

## Adj\_SE3 small adjoint

For `a = (u; ω) ∈ se(3)`, the small adjoint (also called the Lie-algebra
adjoint or the hat map) is:

```
ad(a) = [ω̃   ũ]   ∈ ℝ⁶ˣ⁶
        [0    ω̃]
```

`ad(a) · b = [a, b]` is the Lie bracket of `a` and `b`.  This matrix is
identical to [`sk_SE3`](@ref) applied to `a`.

## Adj\_SE3 Lie bracket

For `a = (u_a; ω_a)` and `b = (u_b; ω_b)`, the Lie bracket computed
matrix-free is:

```
[a, b][1:3] = ω_a × u_b  +  u_a × ω_b
[a, b][4:6] = ω_a × ω_b
```

The fourth overload `Adj_SE3(a, b)` returns this vector directly without
forming the 6×6 matrix.

## Adj\_SE3 methods

| Signature                           | Returns  | Description                        |
|-------------------------------------|----------|------------------------------------|
| `Adj_SE3(H::NodeFrame)`             | `MAT6`   | Full 6×6 group adjoint `Ad(H)`     |
| `Adj_SE3(H::NodeFrame, a::VEC6)`    | `VEC6`   | Matrix-free `Ad(H)·a`              |
| `Adj_SE3(a::VEC6)`                  | `MAT6`   | Small adjoint `ad(a)` = `sk_SE3(a)`|
| `Adj_SE3(a::VEC6, b::VEC6)`         | `VEC6`   | Lie bracket `[a, b]`               |

**Dependencies:** [`R_SO3`](@ref), [`tilde`](@ref), [`crossp`](@ref)

**See also:** [`sk_SE3`](@ref), [`exp_SE3`](@ref), [`NodeFrame`](@ref)

```@docs
Adj_SE3
```
