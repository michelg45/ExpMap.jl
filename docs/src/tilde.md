```@meta
CurrentModule = ExpMap
```

# tilde — Skew-Symmetric (Hat) Map

`tilde(v)` constructs the 3×3 **skew-symmetric matrix** associated with a 3D
vector, also known as the **hat map** or **cross-product matrix**.

## tilde definition

Given **v = (vx, vy, vz)ᵀ**:

```
        ┌  0   -vz   vy ┐
ṽ  =   │  vz   0   -vx │
        └ -vy   vx   0  ┘
```

**Key property**: `ṽ · w = v × w` for any vector w ∈ ℝ³.

## Properties

- **Skew-symmetry**: `tilde(v)ᵀ = −tilde(v)`
- **Rank 2** (for v ≠ 0), one zero eigenvalue along v
- **Null space**: `tilde(v) · v = 0`

## tilde role in SO(3)

`tilde` appears in the **Rodrigues rotation formula** ([`R_SO3`](@ref)):

```
R(ψ) = cos θ · I  +  sin θ · k̃  +  (1 − cos θ) · k⊗k
```

and identifies the Lie algebra so(3) with skew-symmetric matrices.

## tilde methods

| Signature              | Description                                     |
|------------------------|-------------------------------------------------|
| `tilde(a::VEC3) → MAT3`| Skew-symmetric matrix of a `VEC3`               |
| `tilde(a::RV3)  → MAT3`| Same for `RV3` (avoids explicit conversion)     |

**See also:** [`outerp`](@ref), [`R_SO3`](@ref)

```@docs
tilde
```
