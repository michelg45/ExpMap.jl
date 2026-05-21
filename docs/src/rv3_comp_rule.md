```@meta
CurrentModule = ExpMap
```

# RV3 Composition Rule

The two-argument constructor `RV3(a, b)` extends the `RV3` type with the
**CRV composition rule**: it returns the Cartesian Rotation Vector `c` such
that R(c) = R(a) · R(b).

## CRV composition background

In SO(3), rotations compose by matrix multiplication:

```
R(c) = R(a) · R(b)
```

This is implemented by routing through the rotation matrix representation:

```
c = invR_SO3( R_SO3(a) · R_SO3(b) )
```

!!! warning "Non-commutativity"
    Composition of CRVs is **non-commutative**: `RV3(a, b) ≠ RV3(b, a)` in general.

## CRV composition uniqueness

The result is unique as long as ‖c‖ < π.  For ‖c‖ ≥ π the logarithmic map
`invR_SO3` may return an equivalent rotation with a different axis orientation.

## CRV composition example

```julia
phi   = RV3(0.0,  2.0, -1.0)
psi   = RV3(5.0,  2.0,  1.0)
theta = RV3(phi, psi)
# → RV3(-6.605006e-01, 1.995766e+00, -2.754112e-01)
```

**Dependencies:** [`R_SO3`](@ref), [`invR_SO3`](@ref)

```@docs
RV3(a::RV3, b::RV3)
```
