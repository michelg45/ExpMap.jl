"""
    T_functions(x::Float64; der::Int = 0) -> Tuple

Compute the scalar coefficient functions `α(x)` and `β(x)` that appear in the
tangent operator `T_SO3`, together with their derivatives up to the requested order.

The tangent operator for SO(3) can be written as

    T(φ) = I + α(‖φ‖) φ̃/‖φ‖ + β(‖φ‖) (φ̂⊗φ̂ − I)

where `φ̂ = φ/‖φ‖` is the unit rotation axis and `φ̃` denotes the skew-symmetric
matrix associated with `φ`.  The scalar functions are

    α(x) = (cos x − 1) / x
    β(x) = 1 − sin(x) / x

For small `x` (‖φ‖ < 1e-3) Taylor series are used to avoid cancellation errors;
trigonometric expressions are used otherwise.

# Arguments
- `x::Float64`  : norm of the rotation vector ‖φ‖ (non-negative).
- `der::Int = 0`: highest derivative order to return (0, 1, or 2).

# Returns
A tuple whose contents depend on `der`:
- `der == 0`: `(α, β)`
- `der == 1`: `(α, α′, β, β′)`
- `der == 2`: `(α, α′, α″, β, β′, β″)`

# Notes
Higher-order derivatives are needed when computing the directional (Gâteaux)
derivative of `T_SO3` with respect to the rotation vector.
"""
function T_functions(x::Float64; der::Int = 0)

    if x > 1.0e-3
        # ── Trigonometric branch ──────────────────────────────────────────────
        sx, cx = sin(x), cos(x)

        α  = (cx - 1) / x
        β  = 1 - sx / x

        if der > 0
            # First derivatives: α′(x), β′(x)
            dα = -sx / x + (1 - cx) / x^2
            dβ =  sx / x^2 - cx / x
        end

        if der > 1
            # Second derivatives: α″(x), β″(x)
            ddα = (-2 + 2cx + 2x*sx - x^2*cx) / x^3
            ddβ =  cx / x^2 + cx / x^2 + sx / x - 2x * sx / x^4
        end

    else
        # ── Taylor-series branch (avoids cancellation for small x) ───────────
        x2, x3, x4, x5, x6 = x^2, x^3, x^4, x^5, x^6

        α  = -(1//2)*x  + (1//24)*x3   - (1//720)*x5
        β  =  (1//6)*x2 - (1//120)*x4  + (1//5040)*x6

        if der > 0
            dα = -(1//2)    + (1//8)*x2  - (1//144)*x4
            dβ =  (1//3)*x  - (1//30)*x3 + (1//840)*x5
        end

        if der > 1
            ddα = (1//4)*x   - (1//36)*x3  + (1//960)*x5
            ddβ = (1//3)     - (1//10)*x2  + (1//168)*x4
        end
    end

    # ── Return the requested set of values ───────────────────────────────────
    der == 0 && return α, β
    der == 1 && return α, dα, β, dβ
    der == 2 && return α, dα, ddα, β, dβ, ddβ
end
