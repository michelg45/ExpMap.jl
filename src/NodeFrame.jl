# ==============================================================================
# NodeFrame.jl — Nodal Frame Type for Julia  (position + orientation)
# ==============================================================================
#
# Defines an immutable struct `NodeFrame` representing the pose (position and
# orientation) of a point or frame in 3D space.
#
# A NodeFrame stores:
#   • x   ::VEC3  — position of the origin in the reference frame
#   • phi ::RV3   — orientation encoded as a Cartesian Rotation Vector (CRV)
#
# Mathematically, a NodeFrame corresponds to an element of the Special
# Euclidean group SE(3): the semidirect product of ℝ³ (translations) and
# SO(3) (rotations).  The group law and inversion are implemented accordingly.
#
# CONSTRUCTION
#   NodeFrame(x::VEC3, phi::RV3)   → from a position and a rotation vector
#   NodeFrame()                    → identity frame  (origin, zero rotation)
#
# GROUP OPERATIONS  (SE(3) structure)
#   h1 * h2                        → composition: apply h2 in the frame of h1
#   -h                             → inverse frame  (h⁻¹ such that h * (-h) = Id)
#
# UTILITIES
#   copy(h)                        → deep copy of the frame
#
# MATHEMATICAL DETAILS
#   Composition  h1 * h2:
#     x   = h1.x + R(h1.phi) · h2.x        (R_SO3: rotate h2.x by h1.phi)
#     phi = compose(h1.phi, h2.phi)         (RV3 two-argument form: CRV composition)
#
#   Inverse  -h:
#     x   = -R(-h.phi) · h.x               (R_SO3: rotate h.x by the inverse rotation)
#     phi = -h.phi                          (negate the rotation vector)
#
# EXTERNAL DEPENDENCIES
#   R_SO3(phi::RV3, v::VEC3) → VEC3        rotate vector v by the SO(3) matrix R(phi)
#   RV3(phi1::RV3, phi2::RV3) → RV3        compose two CRV rotations
#   VEC3, RV3                               see VEC3.jl, RV3.jl
#
# NOTES
#   • NodeFrame is immutable: fields cannot be modified after construction.
#   • The `*` operator implements the SE(3) group product, NOT component-wise
#     multiplication. It is non-commutative: h1 * h2 ≠ h2 * h1 in general.
#   • The unary `-` operator computes the group inverse (not arithmetic negation):
#     h * (-h) = NodeFrame()  (identity frame).
# ==============================================================================


# ------------------------------------------------------------------------------
# Struct definition
# ------------------------------------------------------------------------------

"""
    NodeFrame

Immutable representation of a nodal frame in 3D space, encoding both the
position and orientation of a point or coordinate frame.

Fields:
- `x   ::VEC3` — position of the frame origin in the reference frame
- `phi ::RV3`  — orientation as a Cartesian Rotation Vector (CRV); the
                  rotation angle is `‖phi‖` (radians) and the axis is `phi/‖phi‖`

A `NodeFrame` is an element of SE(3), the Special Euclidean group. Group
composition and inversion are available via `*` and unary `-`.

# Example
```julia
x   = VEC3(0.0, 1.0, 0.0)
phi = RV3(0.0, 2.0, -1.0)
H   = NodeFrame(x, phi)   # frame at position x, rotated by phi

H2  = H * H               # compose H with itself
Hi  = -H                  # inverse of H
```
"""
struct NodeFrame
    x::VEC3      # position of the frame origin
    phi::RV3     # orientation as a Cartesian Rotation Vector
end


# ------------------------------------------------------------------------------
# Constructors
# ------------------------------------------------------------------------------

"""
    NodeFrame() → NodeFrame

Construct the identity frame: origin at zero, no rotation.
Equivalent to the identity element of SE(3).
"""
NodeFrame() = NodeFrame(VEC3(), RV3())


# ------------------------------------------------------------------------------
# Group operations  (SE(3))
# ------------------------------------------------------------------------------

"""
    h1 * h2  →  NodeFrame

Compose two frames according to the SE(3) group law:

    x   = h1.x + R(h1.phi) · h2.x
    phi = compose(h1.phi, h2.phi)

Geometrically: express the pose `h2` in the coordinate system defined by `h1`.
This operation is **non-commutative**: `h1 * h2 ≠ h2 * h1` in general.

Relies on:
- `R_SO3(phi, v)` — rotate vector `v` by the SO(3) rotation matrix R(phi)
- `RV3(phi1, phi2)` — compose two Cartesian Rotation Vectors
"""
function Base.:*(h1::NodeFrame, h2::NodeFrame)
    x   = h1.x + R_SO3(h1.phi, h2.x)   # translate h2's origin into h1's frame
    phi = RV3(h1.phi, h2.phi)            # compose the two orientations
    return NodeFrame(x, phi)
end

"""
    -h  →  NodeFrame

Return the SE(3) group inverse of `h`, i.e. the frame `h⁻¹` such that:

    h * (-h) = NodeFrame()   (identity frame)

    x   = -R(-h.phi) · h.x
    phi = -h.phi

Relies on:
- `R_SO3(phi, v)` — rotate vector `v` by R(phi); called here with `-h.phi`
"""
function Base.:-(h::NodeFrame)
    x   = -R_SO3(-h.phi, h.x)   # un-rotate then negate the position
    phi = -h.phi                  # reverse the rotation
    return NodeFrame(x, phi)
end


# ------------------------------------------------------------------------------
# Utilities
# ------------------------------------------------------------------------------

"""
    copy(h::NodeFrame) → NodeFrame

Return a deep copy of `h` (useful when storing `NodeFrame` inside mutable containers).
"""
Base.copy(h::NodeFrame) = NodeFrame(copy(h.x), copy(h.phi))
