```@meta
CurrentModule = ExpMap
```

# frame\_on\_line — Orthonormal Frame Aligned on a Line Segment

`frame_on_line` constructs an orthonormal frame whose first axis **n₁** is
aligned along a given line segment, and returns the corresponding
**Cartesian Rotation Vector** (CRV).

## Algorithm

Given two endpoint positions `x_1`, `x_2` ∈ ℝ³:

1. Compute the unit tangent **n₁ = (x₂ − x₁) / ‖x₂ − x₁‖**.
2. Pick a random unit vector **r** such that `|n₁ · r| ≤ 0.5`
   (up to 5 attempts), ensuring **r** is not nearly parallel to **n₁**.
3. Orthogonalise via cross products (Gram-Schmidt):

```
n₃ = normalise(n₁ × r)
n₂ = normalise(n₃ × n₁)
```

4. Assemble the rotation matrix **R = [n₁ | n₂ | n₃]** and extract the CRV
   via the logarithmic map [`invR_SO3`](@ref).

The resulting frame satisfies **R · e₁ = n₁** — the first local axis points
along the segment.

## Non-deterministic output

The axes **n₂** and **n₃** depend on the random draw in step 2.
For reproducible results, seed Julia's RNG before calling:

```julia
using Random
Random.seed!(42)
rv = frame_on_line(x_1, x_2)
```

**Dependencies:** [`invR_SO3`](@ref), [`MAT3`](@ref), `norm2`, `dotp`, `crossp`

**See also:** [`invR_SO3`](@ref), [`R_SO3`](@ref), [`NodeFrame`](@ref)

```@docs
frame_on_line
```
