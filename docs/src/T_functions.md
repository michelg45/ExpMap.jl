```@meta
CurrentModule = ExpMap
```

# T_functions — Scalar Coefficients α(x), β(x)

`T_functions` computes the scalar coefficient functions **α(x)** and **β(x)**
that appear in the tangent operator T(φ) of SO(3), together with their
derivatives up to the requested order.

## T\_functions background

The tangent operator can be written as:

```
T(φ) = I + α(‖φ‖) φ̃/‖φ‖ + β(‖φ‖) (φ̂⊗φ̂ − I)
```

The scalar functions are:

```
α(x) = (cos x − 1) / x
β(x) = 1 − sin(x) / x
```

For small x (‖φ‖ < 1e-3) **Taylor series** are used to avoid cancellation
errors; trigonometric expressions are used otherwise.

## T\_functions return values

| `der` | Returned tuple              |
|-------|-----------------------------|
| `0`   | `(α, β)`                    |
| `1`   | `(α, α′, β, β′)`            |
| `2`   | `(α, α′, α″, β, β′, β″)`   |

Higher-order derivatives are needed when computing the directional (Gâteaux)
derivative of T_SO3 with respect to the rotation vector.

**See also:** [`invT_functions`](@ref), [`invT_SO3`](@ref)

```@docs
T_functions
```
