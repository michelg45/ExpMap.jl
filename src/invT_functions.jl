"""
    invT_functions(x::Float64; der::Int = 0) -> Tuple

Compute the scalar coefficient function `γ(x)` that appears in the inverse
tangent operator `invT_SO3`, together with its derivatives up to the requested order.

The inverse tangent operator for SO(3) can be written as

    T⁻¹(φ) = I + ½ φ̃ + γ(‖φ‖) (φ̂⊗φ̂ − I)

where `φ̂ = φ/‖φ‖` is the unit rotation axis, `φ̃` is the skew-symmetric matrix
associated with `φ`, and

    γ(x) = 1 − (x/2) cot(x/2)
          = 1 − (x/2) cos(x/2) / sin(x/2)

For small `x` (‖φ‖ < 1e-3) Taylor series are used to avoid cancellation errors;
trigonometric expressions are used otherwise.

# Arguments
- `x::Float64`  : norm of the rotation vector ‖φ‖ (non-negative).
- `der::Int = 0`: highest derivative order to return (0, 1, or 2).

# Returns
A tuple whose contents depend on `der`:
- `der == 0`: `(γ,)`         — returned as a scalar (not a 1-tuple)
- `der == 1`: `(γ, γ′)`
- `der == 2`: `(γ, γ′, γ″)`

# Notes
Higher-order derivatives are needed when computing the directional (Gâteaux)
derivative of `invT_SO3` with respect to the rotation vector.
"""
function invT_functions(x::Float64; der::Int = 0)

    if x > 1.0e-3
        # ── Trigonometric branch ──────────────────────────────────────────────
        sh, ch = sin(x/2), cos(x/2)        # half-angle sines/cosines

        γ   = 1 - (x/2) * ch / sh

        if der > 0
            # First derivative γ′(x)
            dγ  = (-(1//2)*sh*ch + (1//4)*x*sh^2 + (1//4)*x*ch^2) / sh^2
        end

        if der > 1
            # Second derivative γ″(x)
            ddγ = ((1//4)*sin(x)*ch + (1//2)*sh^3
                   - (1//4)*x*sh^2*ch - (1//4)*x*ch^3) / sh^3
        end

    else
        # ── Taylor-series branch (avoids cancellation for small x) ───────────
        x2, x3, x4 = x^2, x^3, x^4

        γ   =  (1//12)*x2 + (1//720)*x4

        if der > 0
            dγ  = (1//6)*x  + (1//180)*x3
        end

        if der > 1
            ddγ = (1//6) + (1//60)*x2 + (1//1008)*x4
        end
    end

    # ── Return the requested set of values ───────────────────────────────────
    der == 0 && return γ
    der == 1 && return γ, dγ
    der == 2 && return γ, dγ, ddγ
end
