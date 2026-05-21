```@meta
CurrentModule = ExpMap
```

# invT_functions — Scalar Coefficient γ(x)

`invT_functions` computes the scalar coefficient function **γ(x)** that
appears in the inverse tangent operator [`invT_SO3`](@ref), together with
its derivatives up to the requested order.

## invT\_functions background

```
γ(x) = 1 − (x/2) cot(x/2)
      = 1 − (x/2) cos(x/2) / sin(x/2)
```

For small x (‖φ‖ < 1e-3) **Taylor series** are used to avoid cancellation
errors; trigonometric expressions are used otherwise.

## invT\_functions return values

| `der` | Returned value / tuple  |
|-------|-------------------------|
| `0`   | scalar `γ`              |
| `1`   | `(γ, γ′)`               |
| `2`   | `(γ, γ′, γ″)`           |

Higher-order derivatives are needed when computing the directional (Gâteaux)
derivative of `invT_SO3` with respect to the rotation vector.

**See also:** [`T_functions`](@ref), [`invT_SO3`](@ref)

```@docs
invT_functions
```
