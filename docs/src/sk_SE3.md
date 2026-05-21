```@meta
CurrentModule = ExpMap
```

# sk_SE3 — Hat Map for SE(3)

`sk_SE3` implements the **hat map** for SE(3), mapping a 6-component twist
vector `v = (u; ω) ∈ ℝ⁶` to its 6×6 block-skew matrix representation in
the Lie algebra se(3).

## sk\_SE3 background

The twist vector `v` is partitioned as:

| Components | Symbol | Meaning              |
|------------|--------|----------------------|
| `v[1:3]`   | `u`    | Translational part   |
| `v[4:6]`   | `ω`    | Angular (rotational) part |

The hat map formula is:

```
sk_SE3(v) = [ω̃   ũ]   ∈ ℝ⁶ˣ⁶
            [0    ω̃]
```

where `ω̃ = tilde(ω)` and `ũ = tilde(u)` are 3×3 skew-symmetric matrices.

## sk\_SE3 action on a twist

For any `a = (a_u; a_ω) ∈ ℝ⁶`, the matrix action recovers the se(3) Lie bracket:

```
sk_SE3(v) · a = [ω × a_u + u × a_ω]  =  [v, a]
                [ω × a_ω           ]
```

## sk\_SE3 relation to Adj\_SE3

`sk_SE3(v)` and `Adj_SE3(v::VEC6)` return the same 6×6 matrix.  The two names
emphasise different interpretations: `sk_SE3` stresses the hat map, while
[`Adj_SE3`](@ref) stresses the small-adjoint (Lie-algebra adjoint) role.

**Dependencies:** [`tilde`](@ref)

**See also:** [`Adj_SE3`](@ref), [`bracket`](@ref)

```@docs
sk_SE3
```
