```@meta
CurrentModule = ExpMap
```

# T\_SE3 — SE(3) Tangent Operator

`T_SE3` constructs the 6×6 **SE(3) tangent operator** from either a `VEC6`
parameter vector (direct evaluation) or a precomputed [`T_SE3_data`](@ref)
struct (shared-data path).

## T\_SE3 background

For the parameter vector `p = [u; φ]` with translation `u ∈ ℝ³` and rotation
vector `φ ∈ ℝ³`, the SE(3) tangent operator has the block upper-triangular form:

```
T_SE3(p) = [ T(φ)    dTu − tilde(Tu) · T(φ) ]
           [ 0₃ₓ₃    T(φ)                    ]
```

with `T(φ) = T_SO3(φ)`, `Tu = T(φ)·u`, and `dTu = DT_SO3(φ, u)`.
The off-diagonal block encodes the coupling between translation and rotation.

## T\_SE3 at −p

With `sign_p = "-"` (requires [`T_SE3_data`](@ref) input), the operator is
evaluated at `−p = [−u; −φ]` using the identity **T(−φ) = T(φ)ᵀ**:

```
T(φ)       →   T(φ)ᵀ
Tu         →  −T(φ)ᵀ·u
dTu[:,i]   →   (∂T/∂φᵢ)ᵀ·u
```

## T\_SE3 action on a variation

```
T_SE3(p) · [δu; δφ] = [ T·δu + (dTu − tilde(Tu)·T)·δφ ]
                       [ T·δφ                            ]
```

The lower block is the standard SO(3) tangent map; the upper block mixes
translational and rotational variations.

## T\_SE3 example

```julia
p = VEC6(1.0, 0.0, 0.0,  0.0, 0.0, π/4)
M = T_SE3(p)                    # 6×6 MAT6

# With precomputed data (shares work with DT_SE3, DinvT_SE3)
a  = T_SE3_input_data(p)
M  = T_SE3(a, "+")              # T_SE3(p)
Mm = T_SE3(a, "-")              # T_SE3(−p)
```

**Dependencies:** [`T_SE3_input`](@ref)

**See also:** [`invT_SE3`](@ref), [`T_SE3_input_data`](@ref), [`DT_SE3`](@ref)

```@docs
T_SE3
```
