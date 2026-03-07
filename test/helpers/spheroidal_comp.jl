# Generic backend for spheroidal vector wave functions M = ∇ × (ψ a)
# where a is a unit director vector expressed in the spheroidal orthonormal basis.
#
# Strategy: compute M from its definition rather than closed-form expansions,
# to serve as an independent check against the hardcoded formulas.

# ---------------------------------------------------------------------------
# Metric data
# ---------------------------------------------------------------------------

"""
    metric_data_spheroidal(oblate, d, ξ, η[, ϕ])

Return metric scale factors and their first partial derivatives for prolate
(oblate=false) or oblate (oblate=true) spheroidal coordinates.

Prolate:  f1 = ξ²−1,  f2 = ξ²−η²,  f3 = 1−η²
Oblate:   f1 = ξ²+1,  f2 = ξ²+η²,  f3 = 1−η²

hξ = d√(f2/f1),  hη = d√(f2/f3),  hϕ = d√(f1·f3)

Returns `(h, dh)` where
- h  = (hξ, hη, hϕ)
- dh = (∂ξhξ, ∂ηhξ, ∂ϕhξ, ∂ξhη, ∂ηhη, ∂ϕhη, ∂ξhϕ, ∂ηhϕ, ∂ϕhϕ)
  as a NamedTuple.
"""
function metric_data_spheroidal(oblate::Bool, d, ξ, η, ϕ=zero(ξ))
    ξ2 = ξ^2
    η2 = η^2

    f1 = oblate ? (ξ2 + 1)  : (ξ2 - 1)
    f2 = oblate ? (ξ2 + η2) : (ξ2 - η2)
    f3 = 1 - η2

    # ∂f2/∂η differs by sign between prolate and oblate
    ∂ηf2 = oblate ? 2η : -2η

    hξ = d * sqrt(f2 / f1)
    hη = d * sqrt(f2 / f3)
    hϕ = d * sqrt(f1 * f3)

    # ∂ξ hξ: same formula for both (∂ξf1 = ∂ξf2 = 2ξ cancel into f1−f2)
    ∂ξhξ = d * ξ * (f1 - f2) / (f1^(3/2) * sqrt(f2))

    # ∂η hξ: from ∂η sqrt(f2/f1) = ∂ηf2 / (2 √f1 √f2)  (f1 doesn't depend on η)
    ∂ηhξ = d * ∂ηf2 / (2 * sqrt(f1) * sqrt(f2))

    # ∂ξ hη: from ∂ξ sqrt(f2/f3) = 2ξ / (2 √f2 √f3)  (f3 doesn't depend on ξ)
    ∂ξhη = d * ξ / (sqrt(f3) * sqrt(f2))

    # ∂η hη: d·(∂ηf2/(2√f2·√f3) + η·√f2/f3^(3/2))
    # Simplifies to d·η·f1/(f3^(3/2)·√f2) for both cases (verified algebraically).
    ∂ηhη = d * (∂ηf2 / (2 * sqrt(f2) * sqrt(f3)) + η * sqrt(f2) / f3^(3/2))

    # ∂ξ hϕ: ∂ξ√(f1·f3) = ∂ξf1·f3/(2√(f1·f3)) = 2ξ·√f3/(2√f1)
    ∂ξhϕ = d * ξ * sqrt(f3) / sqrt(f1)

    # ∂η hϕ: ∂η√(f1·f3) = f1·∂ηf3/(2√(f1·f3)) = f1·(−2η)/(2√f1·√f3)
    ∂ηhϕ = -d * η * sqrt(f1) / sqrt(f3)

    # No ϕ dependence (axisymmetric coordinate system)
    h = (hξ=hξ, hη=hη, hϕ=hϕ)
    dh = (
        ∂ξhξ=∂ξhξ, ∂ηhξ=∂ηhξ, ∂ϕhξ=zero(hξ),
        ∂ξhη=∂ξhη, ∂ηhη=∂ηhη, ∂ϕhη=zero(hη),
        ∂ξhϕ=∂ξhϕ, ∂ηhϕ=∂ηhϕ, ∂ϕhϕ=zero(hϕ),
    )
    return h, dh
end


# ---------------------------------------------------------------------------
# Director data: ẑ
# ---------------------------------------------------------------------------

"""
    director_data_z(oblate, d, ξ, η[, ϕ])

Components of ẑ in the spheroidal orthonormal basis and their first
partial derivatives.

Prolate/oblate (using f1,f2,f3 as in `metric_data_spheroidal`):

    aξ = η √(f1/f2),   aη = ξ √(f3/f2),   aϕ = 0

Returns `(a, da)` as NamedTuples.
"""
function director_data_z(oblate::Bool, d, ξ, η, ϕ=zero(ξ))
    ξ2 = ξ^2
    η2 = η^2

    f1 = oblate ? (ξ2 + 1)  : (ξ2 - 1)
    f2 = oblate ? (ξ2 + η2) : (ξ2 - η2)
    f3 = 1 - η2

    # ∂f2/∂η
    ∂ηf2 = oblate ? 2η : -2η

    aξ = η * sqrt(f1 / f2)
    aη = ξ * sqrt(f3 / f2)
    aϕ = zero(aξ)

    # ∂ξ aξ: η · ∂ξ√(f1/f2) = η·(∂ξf1/(2√f1·√f2) − √f1·∂ξf2/(2·f2^(3/2)))
    #        ∂ξf1 = ∂ξf2 = 2ξ  for both cases
    ∂ξaξ = η * (ξ / (sqrt(f1) * sqrt(f2)) - ξ * sqrt(f1) / f2^(3/2))

    # ∂η aξ: √(f1/f2) + η·∂η√(f1/f2)
    #        ∂η√(f1/f2) = −√f1·∂ηf2/(2·f2^(3/2))  (f1 independent of η)
    ∂ηaξ = sqrt(f1 / f2) - η * sqrt(f1) * ∂ηf2 / (2 * f2^(3/2))

    # ∂ξ aη: √(f3/f2) + ξ·∂ξ√(f3/f2)
    #        ∂ξ√(f3/f2) = −√f3·∂ξf2/(2·f2^(3/2)) = −ξ·√f3/f2^(3/2)
    ∂ξaη = sqrt(f3 / f2) - ξ^2 * sqrt(f3) / f2^(3/2)

    # ∂η aη: ξ·∂η√(f3/f2) = ξ·(∂ηf3/(2√f3·√f2) − √f3·∂ηf2/(2·f2^(3/2)))
    #        ∂ηf3 = −2η
    ∂ηaη = ξ * (-η / (sqrt(f3) * sqrt(f2)) - sqrt(f3) * ∂ηf2 / (2 * f2^(3/2)))

    a = (aξ=aξ, aη=aη, aϕ=aϕ)
    da = (
        ∂ξaξ=∂ξaξ, ∂ηaξ=∂ηaξ, ∂ϕaξ=zero(aξ),
        ∂ξaη=∂ξaη, ∂ηaη=∂ηaη, ∂ϕaη=zero(aη),
        ∂ξaϕ=zero(aϕ), ∂ηaϕ=zero(aϕ), ∂ϕaϕ=zero(aϕ),
    )
    return a, da
end


function director_data_x(oblate::Bool, d, ξ, η, ϕ=zero(ξ))
    ξ2 = ξ^2
    η2 = η^2

    f1 = oblate ? (ξ2 + 1)  : (ξ2 - 1)
    f2 = oblate ? (ξ2 + η2) : (ξ2 - η2)
    f3 = 1 - η2
    ∂ηf2 = oblate ? 2η : -2η

    A = sqrt(f1 / f2)   # sqrt((ξ²±1)/(ξ²∓η²))
    B = sqrt(f3 / f2)   # sqrt((1-η²)/(ξ²∓η²))

    ∂ξA = ξ/(sqrt(f1)*sqrt(f2)) - ξ*sqrt(f1)/f2^(3/2)
    ∂ηA = -sqrt(f1)*∂ηf2/(2*f2^(3/2))

    ∂ξB = -ξ*sqrt(f3)/f2^(3/2)
    ∂ηB = -η/(sqrt(f3)*sqrt(f2)) - sqrt(f3)*∂ηf2/(2*f2^(3/2))

    cϕ = cos(ϕ)
    sϕ = sin(ϕ)

    aξ =  ξ * B * cϕ
    aη = -η * A * cϕ
    aϕ = -sϕ

    ∂ξaξ = cϕ * (B + ξ*∂ξB)
    ∂ηaξ = cϕ * (ξ*∂ηB)
    ∂ϕaξ = -ξ * B * sϕ

    ∂ξaη = -η * cϕ * ∂ξA
    ∂ηaη = -cϕ * (A + η*∂ηA)
    ∂ϕaη =  η * A * sϕ

    a = (aξ=aξ, aη=aη, aϕ=aϕ)
    da = (
        ∂ξaξ=∂ξaξ, ∂ηaξ=∂ηaξ, ∂ϕaξ=∂ϕaξ,
        ∂ξaη=∂ξaη, ∂ηaη=∂ηaη, ∂ϕaη=∂ϕaη,
        ∂ξaϕ=zero(aϕ), ∂ηaϕ=zero(aϕ), ∂ϕaϕ=-cϕ,
    )
    return a, da
end

function director_data_y(oblate::Bool, d, ξ, η, ϕ=zero(ξ))
    ξ2 = ξ^2
    η2 = η^2

    f1 = oblate ? (ξ2 + 1)  : (ξ2 - 1)
    f2 = oblate ? (ξ2 + η2) : (ξ2 - η2)
    f3 = 1 - η2
    ∂ηf2 = oblate ? 2η : -2η

    A = sqrt(f1 / f2)
    B = sqrt(f3 / f2)

    ∂ξA = ξ/(sqrt(f1)*sqrt(f2)) - ξ*sqrt(f1)/f2^(3/2)
    ∂ηA = -sqrt(f1)*∂ηf2/(2*f2^(3/2))

    ∂ξB = -ξ*sqrt(f3)/f2^(3/2)
    ∂ηB = -η/(sqrt(f3)*sqrt(f2)) - sqrt(f3)*∂ηf2/(2*f2^(3/2))

    cϕ = cos(ϕ)
    sϕ = sin(ϕ)

    aξ =  ξ * B * sϕ
    aη = -η * A * sϕ
    aϕ =  cϕ

    ∂ξaξ = sϕ * (B + ξ*∂ξB)
    ∂ηaξ = sϕ * (ξ*∂ηB)
    ∂ϕaξ =  ξ * B * cϕ

    ∂ξaη = -η * sϕ * ∂ξA
    ∂ηaη = -sϕ * (A + η*∂ηA)
    ∂ϕaη = -η * A * cϕ

    a = (aξ=aξ, aη=aη, aϕ=aϕ)
    da = (
        ∂ξaξ=∂ξaξ, ∂ηaξ=∂ηaξ, ∂ϕaξ=∂ϕaξ,
        ∂ξaη=∂ξaη, ∂ηaη=∂ηaη, ∂ϕaη=∂ϕaη,
        ∂ξaϕ=zero(aϕ), ∂ηaϕ=zero(aϕ), ∂ϕaϕ=-sϕ,
    )
    return a, da
end

function director_data_r(oblate::Bool, d, ξ, η, ϕ=zero(ξ))
    ξ2 = ξ^2
    η2 = η^2

    f1 = oblate ? (ξ2 + 1)  : (ξ2 - 1)
    f2 = oblate ? (ξ2 + η2) : (ξ2 - η2)
    f3 = 1 - η2
    ∂ηf2 = oblate ? 2η : -2η

    A = sqrt(f1 / f2)
    B = sqrt(f3 / f2)

    ∂ξA = ξ/(sqrt(f1)*sqrt(f2)) - ξ*sqrt(f1)/f2^(3/2)
    ∂ηA = -sqrt(f1)*∂ηf2/(2*f2^(3/2))

    ∂ξB = -ξ*sqrt(f3)/f2^(3/2)
    ∂ηB = -η/(sqrt(f3)*sqrt(f2)) - sqrt(f3)*∂ηf2/(2*f2^(3/2))

    s = oblate ? -1 : +1
    c = 0.5 * d

    aξ = c * ξ * A
    aη = c * s * η * B
    aϕ = zero(aξ)

    ∂ξaξ = c * (A + ξ*∂ξA)
    ∂ηaξ = c * (ξ*∂ηA)

    ∂ξaη = c * s * (η*∂ξB)
    ∂ηaη = c * s * (B + η*∂ηB)

    a = (aξ=aξ, aη=aη, aϕ=aϕ)
    da = (
        ∂ξaξ=∂ξaξ, ∂ηaξ=∂ηaξ, ∂ϕaξ=zero(aξ),
        ∂ξaη=∂ξaη, ∂ηaη=∂ηaη, ∂ϕaη=zero(aη),
        ∂ξaϕ=zero(aϕ), ∂ηaϕ=zero(aϕ), ∂ϕaϕ=zero(aϕ),
    )
    return a, da
end

# ---------------------------------------------------------------------------
# Scalar mode: ψ = S(η) R(ξ) cos/sin(mϕ) and its partial derivatives
# ---------------------------------------------------------------------------

"""
    scalar_mode_first_derivatives(R, dR, S, dS, m, ϕ; even=true)

Build the real-parity scalar ψ and its partial derivatives (∂ξ, ∂η, ∂ϕ) from
the separated amplitudes R, S and their ξ, η derivatives dR, dS.

even=true  →  ψ = S R cos(mϕ)
even=false →  ψ = S R sin(mϕ)
"""
function scalar_mode_first_derivatives(R, dR, S, dS, m, ϕ; even::Bool=true)
    smϕ, cmϕ = sincos(m * ϕ)
    ang  = even ? cmϕ : smϕ           # value of angular factor
    dang = even ? (-m * smϕ) : (m * cmϕ)  # ∂ϕ of angular factor

    ψ   = S  * R  * ang
    dψξ = S  * dR * ang
    dψη = dS * R  * ang
    dψϕ = S  * R  * dang

    return (ψ=ψ, dψξ=dψξ, dψη=dψη, dψϕ=dψϕ)
end


# ---------------------------------------------------------------------------
# Generic curl: M = ∇ × (ψ a)
# ---------------------------------------------------------------------------

"""
    curl_first_from_scalar(ψ, dψξ, dψη, dψϕ, h, dh, a, da)

Compute M = ∇ × (ψ **a**) in an orthogonal curvilinear coordinate system
(ξ, η, ϕ) with scale factors h = (hξ, hη, hϕ).

Uses the identity:
    [∇×(ψa)]_i = (1/(h_j h_k)) * (∂_j(h_k ψ a_k) − ∂_k(h_j ψ a_j))

with all metric and director derivatives supplied through `dh` and `da`.
Valid for any director **a** (not restricted to ẑ).
"""
function curl_first_from_scalar(ψ, dψξ, dψη, dψϕ, h, dh, a, da)
    hξ, hη, hϕ = h.hξ, h.hη, h.hϕ
    aξ, aη, aϕ = a.aξ, a.aη, a.aϕ

    ∂ξhξ, ∂ηhξ, ∂ϕhξ = dh.∂ξhξ, dh.∂ηhξ, dh.∂ϕhξ
    ∂ξhη, ∂ηhη, ∂ϕhη = dh.∂ξhη, dh.∂ηhη, dh.∂ϕhη
    ∂ξhϕ, ∂ηhϕ, ∂ϕhϕ = dh.∂ξhϕ, dh.∂ηhϕ, dh.∂ϕhϕ

    ∂ξaξ, ∂ηaξ, ∂ϕaξ = da.∂ξaξ, da.∂ηaξ, da.∂ϕaξ
    ∂ξaη, ∂ηaη, ∂ϕaη = da.∂ξaη, da.∂ηaη, da.∂ϕaη
    ∂ξaϕ, ∂ηaϕ, ∂ϕaϕ = da.∂ξaϕ, da.∂ηaϕ, da.∂ϕaϕ

    # Mξ = (1/(hη hϕ)) * (∂η(hϕ ψ aϕ) − ∂ϕ(hη ψ aη))
    ∂η_hϕψaϕ = (∂ηhϕ * aϕ + hϕ * ∂ηaϕ) * ψ + hϕ * aϕ * dψη
    ∂ϕ_hηψaη = (∂ϕhη * aη + hη * ∂ϕaη) * ψ + hη * aη * dψϕ
    Mξ = (∂η_hϕψaϕ - ∂ϕ_hηψaη) / (hη * hϕ)

    # Mη = (1/(hϕ hξ)) * (∂ϕ(hξ ψ aξ) − ∂ξ(hϕ ψ aϕ))
    ∂ϕ_hξψaξ = (∂ϕhξ * aξ + hξ * ∂ϕaξ) * ψ + hξ * aξ * dψϕ
    ∂ξ_hϕψaϕ = (∂ξhϕ * aϕ + hϕ * ∂ξaϕ) * ψ + hϕ * aϕ * dψξ
    Mη = (∂ϕ_hξψaξ - ∂ξ_hϕψaϕ) / (hϕ * hξ)

    # Mϕ = (1/(hξ hη)) * (∂ξ(hη ψ aη) − ∂η(hξ ψ aξ))
    ∂ξ_hηψaη = (∂ξhη * aη + hη * ∂ξaη) * ψ + hη * aη * dψξ
    ∂η_hξψaξ = (∂ηhξ * aξ + hξ * ∂ηaξ) * ψ + hξ * aξ * dψη
    Mϕ = (∂ξ_hηψaη - ∂η_hξψaξ) / (hξ * hη)

    return Mξ, Mη, Mϕ
end


# ---------------------------------------------------------------------------
# M^z from the generic backend
# ---------------------------------------------------------------------------

function _director_data_from_symbol(comp::Symbol, oblate::Bool, d, ξ, η, ϕ)
    if comp == :x
        return director_data_x(oblate, d, ξ, η, ϕ)
    elseif comp == :y
        return director_data_y(oblate, d, ξ, η, ϕ)
    elseif comp == :z
        return director_data_z(oblate, d, ξ, η, ϕ)
    elseif comp == :r
        return director_data_r(oblate, d, ξ, η, ϕ)
    else
        throw(ArgumentError("Unsupported component symbol: $comp. Use :x, :y, :z, or :r"))
    end
end

"""
    M_from_mode(comp::Symbol, m, c, ξ, η, ϕ, k, basis; even=true, oblate=false)

Compute the vector spheroidal wave function:
    M = ∇ × (ψ a)
where `a` is selected by `comp`:
- `:x` -> x̂
- `:y` -> ŷ
- `:z` -> ẑ
- `:r` -> r
"""
function M_from_mode(comp::Symbol, m, c, ξ, η, ϕ, k, basis; even::Bool=true, oblate::Bool=false)
    R, dR, _ = evaluate_radial4(basis, ξ)
    S, dS, _ = evaluate_angular(basis, η)

    ψdata = scalar_mode_first_derivatives(R, dR, S, dS, m, ϕ; even=even)

    d = c / k
    h, dh = metric_data_spheroidal(oblate, d, ξ, η, ϕ)
    a, da = _director_data_from_symbol(comp, oblate, d, ξ, η, ϕ)

    return curl_first_from_scalar(
        ψdata.ψ, ψdata.dψξ, ψdata.dψη, ψdata.dψϕ,
        h, dh, a, da,
    )
end

Mz_from_mode(m, c, ξ, η, ϕ, k, basis; even::Bool=true, oblate::Bool=false) =
    M_from_mode(:z, m, c, ξ, η, ϕ, k, basis; even=even, oblate=oblate)
Mx_from_mode(m, c, ξ, η, ϕ, k, basis; even::Bool=true, oblate::Bool=false) =
    M_from_mode(:x, m, c, ξ, η, ϕ, k, basis; even=even, oblate=oblate)
My_from_mode(m, c, ξ, η, ϕ, k, basis; even::Bool=true, oblate::Bool=false) =
    M_from_mode(:y, m, c, ξ, η, ϕ, k, basis; even=even, oblate=oblate)
Mr_from_mode(m, c, ξ, η, ϕ, k, basis; even::Bool=true, oblate::Bool=false) =
    M_from_mode(:r, m, c, ξ, η, ϕ, k, basis; even=even, oblate=oblate)


function _hardcoded_M_from_symbol(comp::Symbol, ξ, η, ϕ, S, dS, R, dR, m, d, even, type)
    if comp == :x
        return Mˣ(ξ, η, ϕ, S, dS, R, dR, m, d, even, type)
    elseif comp == :y
        return Mʸ(ξ, η, ϕ, S, dS, R, dR, m, d, even, type)
    elseif comp == :z
        return Mᶻ(ξ, η, ϕ, S, dS, R, dR, m, d, even, type)
    elseif comp == :r
        return Mʳ(ξ, η, ϕ, S, dS, R, dR, m, d, even, type)
    else
        throw(ArgumentError("Unsupported component symbol: $comp. Use :x, :y, :z, or :r"))
    end
end

# ---------------------------------------------------------------------------
# Comparison helper
# ---------------------------------------------------------------------------

"""
    compare_M_generic_hardcoded(comp::Symbol, m, c, ξ, η, ϕ, k, basis; even=true, oblate=false, atol=1e-12)

Compare generic `M_from_mode(comp, ...)` against hardcoded `M{comp}`.
"""
function compare_M_generic_hardcoded(comp::Symbol, m, c, ξ, η, ϕ, k, basis;
        even::Bool=true, oblate::Bool=false, atol=1e-12)
    R, dR, _ = evaluate_radial4(basis, ξ)
    S, dS, _ = evaluate_angular(basis, η)

    d = c / k
    type = oblate ? :oblate : :prolate

    gMξ, gMη, gMϕ = M_from_mode(comp, m, c, ξ, η, ϕ, k, basis; even=even, oblate=oblate)
    hMξ, hMη, hMϕ = _hardcoded_M_from_symbol(comp, ξ, η, ϕ, S, dS, R, dR, m, d, even, type)

    Δξ = gMξ - hMξ
    Δη = gMη - hMη
    Δϕ = gMϕ - hMϕ

    relξ = abs(Δξ) / max(abs(gMξ), abs(hMξ), atol)
    relη = abs(Δη) / max(abs(gMη), abs(hMη), atol)
    relϕ = abs(Δϕ) / max(abs(gMϕ), abs(hMϕ), atol)

    return (
        generic   = (ξ=gMξ, η=gMη, ϕ=gMϕ),
        hardcoded = (ξ=hMξ, η=hMη, ϕ=hMϕ),
        abs_error = (ξ=abs(Δξ), η=abs(Δη), ϕ=abs(Δϕ)),
        rel_error = (ξ=relξ, η=relη, ϕ=relϕ),
    )
end

compare_Mz_generic_hardcoded(m, c, ξ, η, ϕ, k, basis; even::Bool=true, oblate::Bool=false, atol=1e-12) =
    compare_M_generic_hardcoded(:z, m, c, ξ, η, ϕ, k, basis; even=even, oblate=oblate, atol=atol)
compare_Mx_generic_hardcoded(m, c, ξ, η, ϕ, k, basis; even::Bool=true, oblate::Bool=false, atol=1e-12) =
    compare_M_generic_hardcoded(:x, m, c, ξ, η, ϕ, k, basis; even=even, oblate=oblate, atol=atol)
compare_My_generic_hardcoded(m, c, ξ, η, ϕ, k, basis; even::Bool=true, oblate::Bool=false, atol=1e-12) =
    compare_M_generic_hardcoded(:y, m, c, ξ, η, ϕ, k, basis; even=even, oblate=oblate, atol=atol)
compare_Mr_generic_hardcoded(m, c, ξ, η, ϕ, k, basis; even::Bool=true, oblate::Bool=false, atol=1e-12) =
    compare_M_generic_hardcoded(:r, m, c, ξ, η, ϕ, k, basis; even=even, oblate=oblate, atol=atol)



fit_alpha(h, g) = dot(g, h) / dot(g, g)

cosine_sim(h, g) = abs(dot(h, g)) / (norm(h) * norm(g) + eps(Float64))

shape_err(h, g) = begin
      α = fit_alpha(h, g)
      norm(h .- α .* g) / (norm(h) + eps(Float64))
  end

function compare_M_alpha_corr_sweep(comp::Symbol, m, c, η, ϕ, k, basis;
        even::Bool=true, oblate::Bool=false, ξs=nothing)

    ξs === nothing && (ξs = oblate ? range(0.05, 3.0, length=200) : range(1.05, 3.0, length=200))

    hξ = ComplexF64[]; hη = ComplexF64[]; hϕ = ComplexF64[]
    gξ = ComplexF64[]; gη = ComplexF64[]; gϕ = ComplexF64[]

    for ξ in ξs
        r = compare_M_generic_hardcoded(comp, m, c, ξ, η, ϕ, k, basis; even=even, oblate=oblate)
        push!(hξ, r.hardcoded.ξ); push!(hη, r.hardcoded.η); push!(hϕ, r.hardcoded.ϕ)
        push!(gξ, r.generic.ξ);   push!(gη, r.generic.η);   push!(gϕ, r.generic.ϕ)
    end

    return (
        ξ = (alpha = fit_alpha(hξ, gξ), corr = cosine_sim(hξ, gξ), err = shape_err(hξ, gξ)),
        η = (alpha = fit_alpha(hη, gη), corr = cosine_sim(hη, gη), err = shape_err(hη, gη)),
        ϕ = (alpha = fit_alpha(hϕ, gϕ), corr = cosine_sim(hϕ, gϕ), err = shape_err(hϕ, gϕ)),
    )
end

compare_Mz_alpha_corr_sweep(m, c, η, ϕ, k, basis; even::Bool=true, oblate::Bool=false, ξs=nothing) =
    compare_M_alpha_corr_sweep(:z, m, c, η, ϕ, k, basis; even=even, oblate=oblate, ξs=ξs)
compare_Mx_alpha_corr_sweep(m, c, η, ϕ, k, basis; even::Bool=true, oblate::Bool=false, ξs=nothing) =
    compare_M_alpha_corr_sweep(:x, m, c, η, ϕ, k, basis; even=even, oblate=oblate, ξs=ξs)
compare_My_alpha_corr_sweep(m, c, η, ϕ, k, basis; even::Bool=true, oblate::Bool=false, ξs=nothing) =
    compare_M_alpha_corr_sweep(:y, m, c, η, ϕ, k, basis; even=even, oblate=oblate, ξs=ξs)
compare_Mr_alpha_corr_sweep(m, c, η, ϕ, k, basis; even::Bool=true, oblate::Bool=false, ξs=nothing) =
    compare_M_alpha_corr_sweep(:r, m, c, η, ϕ, k, basis; even=even, oblate=oblate, ξs=ξs)



