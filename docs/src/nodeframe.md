```@meta
CurrentModule = ExpMap
```

# NodeFrame — Nodal Frame (SE(3))

`NodeFrame` is an **immutable** struct representing the **pose** (position and
orientation) of a point or coordinate frame in 3D space.

A `NodeFrame` is an element of **SE(3)**, the Special Euclidean group —
the semidirect product of ℝ³ (translations) and SO(3) (rotations).
The group law and inversion are implemented accordingly.

## Fields

| Field  | Type    | Description                                              |
|--------|---------|----------------------------------------------------------|
| `x`    | `VEC3`  | Position of the frame origin in the reference frame       |
| `phi`  | `RV3`   | Orientation as a CRV; ‖phi‖ is the angle, phi/‖phi‖ the axis |

## NodeFrame construction

| Syntax                      | Description                                     |
|-----------------------------|-------------------------------------------------|
| `NodeFrame(x::VEC3, phi::RV3)` | From a position and a rotation vector         |
| `NodeFrame()`               | Identity frame (origin, zero rotation)          |

## Group operations — SE(3)

!!! warning "Non-commutativity"
    `*` implements the **SE(3) group product**, not component-wise multiplication.
    It is **non-commutative**: `h1 * h2 ≠ h2 * h1` in general.

### Composition  `h1 * h2`

Expresses the pose `h2` in the coordinate system defined by `h1`:

```
x   = h1.x + R(h1.phi) · h2.x
phi = RV3(h1.phi, h2.phi)          # CRV composition rule
```

Relies on [`R_SO3(phi, v)`](@ref) (rotate vector `v` by R(phi)) and
[`RV3(a, b)`](@ref) (compose two CRVs).

### Inverse  `-h`

Returns the frame `h⁻¹` such that `h * (-h) = NodeFrame()`:

```
x   = -R(-h.phi) · h.x
phi = -h.phi
```

!!! note
    Unary `-` computes the **group inverse**, not arithmetic negation.

## Utilities

| Function    | Description   |
|-------------|---------------|
| `copy(h)`   | Deep copy      |

## NodeFrame example

```julia
x   = VEC3(0.0, 1.0, 0.0)
phi = RV3(0.0, 2.0, -1.0)
H   = NodeFrame(x, phi)   # frame at position x, rotated by phi

H2  = H * H               # compose H with itself (SE(3) product)
Hi  = -H                  # inverse of H:  H * Hi == NodeFrame()
```

```@docs
NodeFrame
```
