# Structured backend for spheroidal vector wave functions:
#   A = a * ψ
#   M = curl(A)
#   N = (1/k) curl(M)
#
# The implementation is block-based and keeps first/second derivatives together
# to avoid finite differences and external AD packages.

using LinearAlgebra

# ---------------------------------------------------------------------------
# Second-order scalar container in (ξ,η,ϕ)
# ---------------------------------------------------------------------------

struct Diff2{T}
    v::T
    dξ::T
    dη::T
    dϕ::T
    dξξ::T
    dξη::T
    dξϕ::T
    dηη::T
    dηϕ::T
    dϕϕ::T
end

_zero_like(x) = zero(x)
_one_like(x) = one(x)

function _const2(x)
    z = _zero_like(x)
    return Diff2(x, z, z, z, z, z, z, z, z, z)
end

function _var2_ξ(ξ)
    z = _zero_like(ξ)
    o = _one_like(ξ)
    return Diff2(ξ, o, z, z, z, z, z, z, z, z)
end

function _var2_η(η)
    z = _zero_like(η)
    o = _one_like(η)
    return Diff2(η, z, o, z, z, z, z, z, z, z)
end

function _var2_ϕ(ϕ)
    z = _zero_like(ϕ)
    o = _one_like(ϕ)
    return Diff2(ϕ, z, z, o, z, z, z, z, z, z)
end

function _add2(a::Diff2, b::Diff2)
    return Diff2(
        a.v + b.v,
        a.dξ + b.dξ, a.dη + b.dη, a.dϕ + b.dϕ,
        a.dξξ + b.dξξ, a.dξη + b.dξη, a.dξϕ + b.dξϕ,
        a.dηη + b.dηη, a.dηϕ + b.dηϕ, a.dϕϕ + b.dϕϕ,
    )
end

function _sub2(a::Diff2, b::Diff2)
    return Diff2(
        a.v - b.v,
        a.dξ - b.dξ, a.dη - b.dη, a.dϕ - b.dϕ,
        a.dξξ - b.dξξ, a.dξη - b.dξη, a.dξϕ - b.dξϕ,
        a.dηη - b.dηη, a.dηϕ - b.dηϕ, a.dϕϕ - b.dϕϕ,
    )
end

function _scale2(c, a::Diff2)
    return Diff2(
        c * a.v,
        c * a.dξ, c * a.dη, c * a.dϕ,
        c * a.dξξ, c * a.dξη, c * a.dξϕ,
        c * a.dηη, c * a.dηϕ, c * a.dϕϕ,
    )
end

function _mul2(a::Diff2, b::Diff2)
    return Diff2(
        a.v * b.v,
        a.dξ * b.v + a.v * b.dξ,
        a.dη * b.v + a.v * b.dη,
        a.dϕ * b.v + a.v * b.dϕ,
        a.dξξ * b.v + 2 * a.dξ * b.dξ + a.v * b.dξξ,
        a.dξη * b.v + a.dξ * b.dη + a.dη * b.dξ + a.v * b.dξη,
        a.dξϕ * b.v + a.dξ * b.dϕ + a.dϕ * b.dξ + a.v * b.dξϕ,
        a.dηη * b.v + 2 * a.dη * b.dη + a.v * b.dηη,
        a.dηϕ * b.v + a.dη * b.dϕ + a.dϕ * b.dη + a.v * b.dηϕ,
        a.dϕϕ * b.v + 2 * a.dϕ * b.dϕ + a.v * b.dϕϕ,
    )
end

function _inv2(b::Diff2)
    invv = inv(b.v)
    invv2 = invv^2
    invv3 = invv^3

    dξ = -b.dξ * invv2
    dη = -b.dη * invv2
    dϕ = -b.dϕ * invv2

    dξξ = 2 * b.dξ^2 * invv3 - b.dξξ * invv2
    dξη = 2 * b.dξ * b.dη * invv3 - b.dξη * invv2
    dξϕ = 2 * b.dξ * b.dϕ * invv3 - b.dξϕ * invv2
    dηη = 2 * b.dη^2 * invv3 - b.dηη * invv2
    dηϕ = 2 * b.dη * b.dϕ * invv3 - b.dηϕ * invv2
    dϕϕ = 2 * b.dϕ^2 * invv3 - b.dϕϕ * invv2

    return Diff2(invv, dξ, dη, dϕ, dξξ, dξη, dξϕ, dηη, dηϕ, dϕϕ)
end

_div2(a::Diff2, b::Diff2) = _mul2(a, _inv2(b))

function _sqrt2(u::Diff2)
    s = sqrt(u.v)
    s3 = s^3
    c1 = 0.5 / s
    c2 = -0.25 / s3

    dξ = c1 * u.dξ
    dη = c1 * u.dη
    dϕ = c1 * u.dϕ

    dξξ = c1 * u.dξξ + c2 * (u.dξ * u.dξ)
    dξη = c1 * u.dξη + c2 * (u.dξ * u.dη)
    dξϕ = c1 * u.dξϕ + c2 * (u.dξ * u.dϕ)
    dηη = c1 * u.dηη + c2 * (u.dη * u.dη)
    dηϕ = c1 * u.dηϕ + c2 * (u.dη * u.dϕ)
    dϕϕ = c1 * u.dϕϕ + c2 * (u.dϕ * u.dϕ)

    return Diff2(s, dξ, dη, dϕ, dξξ, dξη, dξϕ, dηη, dηϕ, dϕϕ)
end

function _sin2(u::Diff2)
    s = sin(u.v)
    c = cos(u.v)

    dξ = c * u.dξ
    dη = c * u.dη
    dϕ = c * u.dϕ

    dξξ = c * u.dξξ - s * u.dξ^2
    dξη = c * u.dξη - s * u.dξ * u.dη
    dξϕ = c * u.dξϕ - s * u.dξ * u.dϕ
    dηη = c * u.dηη - s * u.dη^2
    dηϕ = c * u.dηϕ - s * u.dη * u.dϕ
    dϕϕ = c * u.dϕϕ - s * u.dϕ^2

    return Diff2(s, dξ, dη, dϕ, dξξ, dξη, dξϕ, dηη, dηϕ, dϕϕ)
end

function _cos2(u::Diff2)
    c = cos(u.v)
    s = sin(u.v)

    dξ = -s * u.dξ
    dη = -s * u.dη
    dϕ = -s * u.dϕ

    dξξ = -s * u.dξξ - c * u.dξ^2
    dξη = -s * u.dξη - c * u.dξ * u.dη
    dξϕ = -s * u.dξϕ - c * u.dξ * u.dϕ
    dηη = -s * u.dηη - c * u.dη^2
    dηϕ = -s * u.dηϕ - c * u.dη * u.dϕ
    dϕϕ = -s * u.dϕϕ - c * u.dϕ^2

    return Diff2(c, dξ, dη, dϕ, dξξ, dξη, dξϕ, dηη, dηϕ, dϕϕ)
end

function _pow2_int(u::Diff2, n::Int)
    n == 0 && return _const2(one(u.v))
    n == 1 && return u
    out = u
    for _ in 2:n
        out = _mul2(out, u)
    end
    return out
end

# ---------------------------------------------------------------------------
# Data structs
# ---------------------------------------------------------------------------

struct Scalar2ndData{T}
    ψ::Diff2{T}
end

struct Metric2ndData{T}
    hξ::Diff2{T}
    hη::Diff2{T}
    hϕ::Diff2{T}
end

struct Director2ndData{T}
    aξ::Diff2{T}
    aη::Diff2{T}
    aϕ::Diff2{T}
end

struct A2ndData{T}
    Aξ::Diff2{T}
    Aη::Diff2{T}
    Aϕ::Diff2{T}
end

struct FirstCurlData
    q
    D
    M
    dq
    dD
    dM
end

# ---------------------------------------------------------------------------
# Scalar ψ data
# ---------------------------------------------------------------------------

function scalar_mode_second_derivatives(R, dR, d2R, S, dS, d2S, m, ϕ; even::Bool=true)
    sm, cm = sincos(m * ϕ)
    ang = even ? cm : sm
    dang = even ? (-m * sm) : (m * cm)
    ddang = even ? (-(m^2) * cm) : (-(m^2) * sm)

    ψ = S * R * ang
    dψξ = S * dR * ang
    dψη = dS * R * ang
    dψϕ = S * R * dang

    dψξξ = S * d2R * ang
    dψηη = d2S * R * ang
    dψϕϕ = S * R * ddang

    dψξη = dS * dR * ang
    dψξϕ = S * dR * dang
    dψηϕ = dS * R * dang

    ψ2 = Diff2(ψ, dψξ, dψη, dψϕ, dψξξ, dψξη, dψξϕ, dψηη, dψηϕ, dψϕϕ)
    return Scalar2ndData(ψ2)
end

# ---------------------------------------------------------------------------
# Metric + director second-order data (analytic via block composition)
# ---------------------------------------------------------------------------

function metric_data_spheroidal_2nd(oblate::Bool, d, ξ, η, ϕ=zero(ξ))
    ξd = _var2_ξ(ξ)
    ηd = _var2_η(η)
    _ = _var2_ϕ(ϕ) # kept for signature symmetry

    one2 = _const2(one(ξ))
    d2 = _const2(d)

    ξ2 = _pow2_int(ξd, 2)
    η2 = _pow2_int(ηd, 2)

    f1 = oblate ? _add2(ξ2, one2) : _sub2(ξ2, one2)
    f2 = oblate ? _add2(ξ2, η2) : _sub2(ξ2, η2)
    f3 = _sub2(one2, η2)

    hξ = _mul2(d2, _sqrt2(_div2(f2, f1)))
    hη = _mul2(d2, _sqrt2(_div2(f2, f3)))
    hϕ = _mul2(d2, _sqrt2(_mul2(f1, f3)))

    return Metric2ndData(hξ, hη, hϕ)
end

function _director_data_symbol_2nd(comp::Symbol, oblate::Bool, d, ξ, η, ϕ=zero(ξ))
    ξd = _var2_ξ(ξ)
    ηd = _var2_η(η)
    ϕd = _var2_ϕ(ϕ)

    one2 = _const2(one(ξ))
    halfd = _const2(0.5 * d)

    ξ2 = _pow2_int(ξd, 2)
    η2 = _pow2_int(ηd, 2)

    f1 = oblate ? _add2(ξ2, one2) : _sub2(ξ2, one2)
    f2 = oblate ? _add2(ξ2, η2) : _sub2(ξ2, η2)
    f3 = _sub2(one2, η2)

    A = _sqrt2(_div2(f1, f2))
    B = _sqrt2(_div2(f3, f2))
    cϕ = _cos2(ϕd)
    sϕ = _sin2(ϕd)

    if comp == :x
        aξ = _mul2(_mul2(ξd, B), cϕ)
        aη = _scale2(-1, _mul2(_mul2(ηd, A), cϕ))
        aϕ = _scale2(-1, sϕ)
    elseif comp == :y
        aξ = _mul2(_mul2(ξd, B), sϕ)
        aη = _scale2(-1, _mul2(_mul2(ηd, A), sϕ))
        aϕ = cϕ
    elseif comp == :z
        aξ = _mul2(ηd, A)
        aη = _mul2(ξd, B)
        aϕ = _const2(zero(ξ))
    elseif comp == :r
        s = oblate ? -1 : 1
        aξ = _mul2(halfd, _mul2(ξd, A))
        aη = _scale2(s, _mul2(halfd, _mul2(ηd, B)))
        aϕ = _const2(zero(ξ))
    else
        throw(ArgumentError("Unsupported component symbol: $comp. Use :x, :y, :z, or :r"))
    end

    return Director2ndData(aξ, aη, aϕ)
end

# ---------------------------------------------------------------------------
# A blocks
# ---------------------------------------------------------------------------

function assemble_A(s::Scalar2ndData, a::Director2ndData)
    return A2ndData(
        _mul2(a.aξ, s.ψ),
        _mul2(a.aη, s.ψ),
        _mul2(a.aϕ, s.ψ),
    )
end

function assemble_dA(A::A2ndData)
    return (
        ξ = (dξ=A.Aξ.dξ, dη=A.Aξ.dη, dϕ=A.Aξ.dϕ),
        η = (dξ=A.Aη.dξ, dη=A.Aη.dη, dϕ=A.Aη.dϕ),
        ϕ = (dξ=A.Aϕ.dξ, dη=A.Aϕ.dη, dϕ=A.Aϕ.dϕ),
    )
end

function assemble_ddA(A::A2ndData)
    return (
        ξ = (dξξ=A.Aξ.dξξ, dξη=A.Aξ.dξη, dξϕ=A.Aξ.dξϕ, dηη=A.Aξ.dηη, dηϕ=A.Aξ.dηϕ, dϕϕ=A.Aξ.dϕϕ),
        η = (dξξ=A.Aη.dξξ, dξη=A.Aη.dξη, dξϕ=A.Aη.dξϕ, dηη=A.Aη.dηη, dηϕ=A.Aη.dηϕ, dϕϕ=A.Aη.dϕϕ),
        ϕ = (dξξ=A.Aϕ.dξξ, dξη=A.Aϕ.dξη, dξϕ=A.Aϕ.dξϕ, dηη=A.Aϕ.dηη, dηϕ=A.Aϕ.dηϕ, dϕϕ=A.Aϕ.dϕϕ),
    )
end

# ---------------------------------------------------------------------------
# First curl, dM, second curl
# ---------------------------------------------------------------------------

_dfirst(x::Diff2) = (dξ=x.dξ, dη=x.dη, dϕ=x.dϕ)

function _q_derivs(HϕAϕ::Diff2, HηAη::Diff2, HξAξ::Diff2)
    # qξ = ∂η(HϕAϕ) - ∂ϕ(HηAη)
    qξ = (
        v  = HϕAϕ.dη  - HηAη.dϕ,
        dξ = HϕAϕ.dξη - HηAη.dξϕ,
        dη = HϕAϕ.dηη - HηAη.dηϕ,
        dϕ = HϕAϕ.dηϕ - HηAη.dϕϕ,
    )
    # qη = ∂ϕ(HξAξ) - ∂ξ(HϕAϕ)
    qη = (
        v  = HξAξ.dϕ  - HϕAϕ.dξ,
        dξ = HξAξ.dξϕ - HϕAϕ.dξξ,
        dη = HξAξ.dηϕ - HϕAϕ.dξη,
        dϕ = HξAξ.dϕϕ - HϕAϕ.dξϕ,
    )
    # qϕ = ∂ξ(HηAη) - ∂η(HξAξ)
    qϕ = (
        v  = HηAη.dξ  - HξAξ.dη,
        dξ = HηAη.dξξ - HξAξ.dξη,
        dη = HηAη.dξη - HξAξ.dηη,
        dϕ = HηAη.dξϕ - HξAξ.dηϕ,
    )
    return (ξ=qξ, η=qη, ϕ=qϕ)
end

_quot_d(qv, dqv, Dv, dDv) = (dqv * Dv - qv * dDv) / (Dv^2)

function first_curl(h::Metric2ndData, A::A2ndData)
    hξ, hη, hϕ = h.hξ, h.hη, h.hϕ
    Aξ, Aη, Aϕ = A.Aξ, A.Aη, A.Aϕ

    # Explicitly use derivatives of blocks for clarity.
    HϕAϕ = _mul2(hϕ, Aϕ)
    HηAη = _mul2(hη, Aη)
    HξAξ = _mul2(hξ, Aξ)

    q = _q_derivs(HϕAϕ, HηAη, HξAξ)

    Dξ = _mul2(hη, hϕ)
    Dη = _mul2(hϕ, hξ)
    Dϕ = _mul2(hξ, hη)
    D = (ξ=Dξ.v, η=Dη.v, ϕ=Dϕ.v)
    dD = (ξ=_dfirst(Dξ), η=_dfirst(Dη), ϕ=_dfirst(Dϕ))

    M = (ξ=q.ξ.v / D.ξ, η=q.η.v / D.η, ϕ=q.ϕ.v / D.ϕ)
    dM = (
        ξ = (
            dξ=_quot_d(q.ξ.v, q.ξ.dξ, D.ξ, dD.ξ.dξ),
            dη=_quot_d(q.ξ.v, q.ξ.dη, D.ξ, dD.ξ.dη),
            dϕ=_quot_d(q.ξ.v, q.ξ.dϕ, D.ξ, dD.ξ.dϕ),
        ),
        η = (
            dξ=_quot_d(q.η.v, q.η.dξ, D.η, dD.η.dξ),
            dη=_quot_d(q.η.v, q.η.dη, D.η, dD.η.dη),
            dϕ=_quot_d(q.η.v, q.η.dϕ, D.η, dD.η.dϕ),
        ),
        ϕ = (
            dξ=_quot_d(q.ϕ.v, q.ϕ.dξ, D.ϕ, dD.ϕ.dξ),
            dη=_quot_d(q.ϕ.v, q.ϕ.dη, D.ϕ, dD.ϕ.dη),
            dϕ=_quot_d(q.ϕ.v, q.ϕ.dϕ, D.ϕ, dD.ϕ.dϕ),
        ),
    )

    dq = (
        ξ = (dξ=q.ξ.dξ, dη=q.ξ.dη, dϕ=q.ξ.dϕ),
        η = (dξ=q.η.dξ, dη=q.η.dη, dϕ=q.η.dϕ),
        ϕ = (dξ=q.ϕ.dξ, dη=q.ϕ.dη, dϕ=q.ϕ.dϕ),
    )
    return FirstCurlData(q, D, M, dq, dD, dM)
end

function differentiate_first_curl(M::FirstCurlData)
    return M.dM
end

function second_curl(h::Metric2ndData, M::FirstCurlData, k)
    hξ, hη, hϕ = h.hξ, h.hη, h.hϕ
    Mξ, Mη, Mϕ = M.M.ξ, M.M.η, M.M.ϕ
    dMξ, dMη, dMϕ = M.dM.ξ, M.dM.η, M.dM.ϕ

    dη_hϕMϕ = hϕ.dη * Mϕ + hϕ.v * dMϕ.dη
    dϕ_hηMη = hη.dϕ * Mη + hη.v * dMη.dϕ
    Cξ = (dη_hϕMϕ - dϕ_hηMη) / (hη.v * hϕ.v)

    dϕ_hξMξ = hξ.dϕ * Mξ + hξ.v * dMξ.dϕ
    dξ_hϕMϕ = hϕ.dξ * Mϕ + hϕ.v * dMϕ.dξ
    Cη = (dϕ_hξMξ - dξ_hϕMϕ) / (hϕ.v * hξ.v)

    dξ_hηMη = hη.dξ * Mη + hη.v * dMη.dξ
    dη_hξMξ = hξ.dη * Mξ + hξ.v * dMξ.dη
    Cϕ = (dξ_hηMη - dη_hξMξ) / (hξ.v * hη.v)

    return (Cξ / k, Cη / k, Cϕ / k)
end

# ---------------------------------------------------------------------------
# High-level API
# ---------------------------------------------------------------------------

"""
    N_from_mode(comp::Symbol, m, c, ξ, η, ϕ, k, basis; even=true, oblate=false)

Block-based evaluation of
    N = (1/k) curl(curl(a ψ))
without finite differences.

`comp` selects the director:
- `:x`, `:y`, `:z`, `:r`
"""
function N_from_mode(comp::Symbol, m, c, ξ, η, ϕ, k, basis; even::Bool=true, oblate::Bool=false)
    R, dR, d2R = evaluate_radial4(basis, ξ)
    S, dS, d2S = evaluate_angular(basis, η)

    sdata = scalar_mode_second_derivatives(R, dR, d2R, S, dS, d2S, m, ϕ; even=even)

    d = c / k
    hdata = metric_data_spheroidal_2nd(oblate, d, ξ, η, ϕ)
    adata = _director_data_symbol_2nd(comp, oblate, d, ξ, η, ϕ)

    A = assemble_A(sdata, adata)
    _ = assemble_dA(A)
    _ = assemble_ddA(A)

    M = first_curl(hdata, A)
    _ = differentiate_first_curl(M)

    return second_curl(hdata, M, k)
end

Nz_from_mode(m, c, ξ, η, ϕ, k, basis; even::Bool=true, oblate::Bool=false) =
    N_from_mode(:z, m, c, ξ, η, ϕ, k, basis; even=even, oblate=oblate)
Nx_from_mode(m, c, ξ, η, ϕ, k, basis; even::Bool=true, oblate::Bool=false) =
    N_from_mode(:x, m, c, ξ, η, ϕ, k, basis; even=even, oblate=oblate)
Ny_from_mode(m, c, ξ, η, ϕ, k, basis; even::Bool=true, oblate::Bool=false) =
    N_from_mode(:y, m, c, ξ, η, ϕ, k, basis; even=even, oblate=oblate)
Nr_from_mode(m, c, ξ, η, ϕ, k, basis; even::Bool=true, oblate::Bool=false) =
    N_from_mode(:r, m, c, ξ, η, ϕ, k, basis; even=even, oblate=oblate)

# ---------------------------------------------------------------------------
# M from block backend + comparisons
# ---------------------------------------------------------------------------

function M_from_mode_blocks(comp::Symbol, m, c, ξ, η, ϕ, k, basis; even::Bool=true, oblate::Bool=false)
    R, dR, d2R = evaluate_radial4(basis, ξ)
    S, dS, d2S = evaluate_angular(basis, η)

    sdata = scalar_mode_second_derivatives(R, dR, d2R, S, dS, d2S, m, ϕ; even=even)
    d = c / k
    hdata = metric_data_spheroidal_2nd(oblate, d, ξ, η, ϕ)
    adata = _director_data_symbol_2nd(comp, oblate, d, ξ, η, ϕ)
    A = assemble_A(sdata, adata)
    M = first_curl(hdata, A)
    return M.M.ξ, M.M.η, M.M.ϕ
end

function _hardcoded_M_from_symbol_local(comp::Symbol, ξ, η, ϕ, S, dS, R, dR, m, d, even, type)
    if comp == :x
        return Mx(ξ, η, ϕ, S, dS, R, dR, m, d, even, type)
    elseif comp == :y
        return My(ξ, η, ϕ, S, dS, R, dR, m, d, even, type)
    elseif comp == :z
        return Mz(ξ, η, ϕ, S, dS, R, dR, m, d, even, type)
    elseif comp == :r
        return Mr(ξ, η, ϕ, S, dS, R, dR, m, d, even, type)
    end
    throw(ArgumentError("Unsupported component symbol: $comp. Use :x, :y, :z, or :r"))
end

function compare_M_blocks_hardcoded(comp::Symbol, m, c, ξ, η, ϕ, k, basis;
        even::Bool=true, oblate::Bool=false, atol=1e-12)
    R, dR, _ = evaluate_radial4(basis, ξ)
    S, dS, _ = evaluate_angular(basis, η)
    d = c / k
    type = oblate ? :oblate : :prolate

    gMξ, gMη, gMϕ = M_from_mode_blocks(comp, m, c, ξ, η, ϕ, k, basis; even=even, oblate=oblate)
    hMξ, hMη, hMϕ = _hardcoded_M_from_symbol_local(comp, ξ, η, ϕ, S, dS, R, dR, m, d, even, type)

    Δξ = gMξ - hMξ
    Δη = gMη - hMη
    Δϕ = gMϕ - hMϕ

    relξ = abs(Δξ) / max(abs(gMξ), abs(hMξ), atol)
    relη = abs(Δη) / max(abs(gMη), abs(hMη), atol)
    relϕ = abs(Δϕ) / max(abs(gMϕ), abs(hMϕ), atol)

    return (
        blocks = (ξ=gMξ, η=gMη, ϕ=gMϕ),
        hardcoded = (ξ=hMξ, η=hMη, ϕ=hMϕ),
        abs_error = (ξ=abs(Δξ), η=abs(Δη), ϕ=abs(Δϕ)),
        rel_error = (ξ=relξ, η=relη, ϕ=relϕ),
    )
end

function compare_M_blocks_alpha_corr_sweep(comp::Symbol, m, c, η, ϕ, k, basis;
        even::Bool=true, oblate::Bool=false, ξs=nothing)
    ξs === nothing && (ξs = oblate ? range(0.05, 3.0, length=200) : range(1.05, 3.0, length=200))

    hξ = ComplexF64[]; hη = ComplexF64[]; hϕ = ComplexF64[]
    gξ = ComplexF64[]; gη = ComplexF64[]; gϕ = ComplexF64[]

    for ξ in ξs
        r = compare_M_blocks_hardcoded(comp, m, c, ξ, η, ϕ, k, basis; even=even, oblate=oblate)
        push!(hξ, r.hardcoded.ξ); push!(hη, r.hardcoded.η); push!(hϕ, r.hardcoded.ϕ)
        push!(gξ, r.blocks.ξ);    push!(gη, r.blocks.η);    push!(gϕ, r.blocks.ϕ)
    end

    return (
        ξ = (alpha = fit_alpha(hξ, gξ), corr = cosine_sim(hξ, gξ), err = shape_err(hξ, gξ)),
        η = (alpha = fit_alpha(hη, gη), corr = cosine_sim(hη, gη), err = shape_err(hη, gη)),
        ϕ = (alpha = fit_alpha(hϕ, gϕ), corr = cosine_sim(hϕ, gϕ), err = shape_err(hϕ, gϕ)),
    )
end

# ---------------------------------------------------------------------------
# Comparisons vs hardcoded references
# ---------------------------------------------------------------------------

fit_alpha(h, g) = dot(g, h) / dot(g, g)
cosine_sim(h, g) = abs(dot(h, g)) / (norm(h) * norm(g) + eps(Float64))
shape_err(h, g) = begin
    α = fit_alpha(h, g)
    norm(h .- α .* g) / (norm(h) + eps(Float64))
end

function _hardcoded_N_from_symbol(comp::Symbol, m, c, ξ, η, ϕ, k, basis, even, oblate)
    if comp == :z
        return Nᶻ(m, c, ξ, η, ϕ, k, basis, even, oblate)
    elseif comp == :x
        return Nˣ(m, c, ξ, η, ϕ, k, basis, even, oblate)
    elseif comp == :y
        return Nʸ(m, c, ξ, η, ϕ, k, basis, even, oblate)
    elseif comp == :r
        return Nʳ(m, c, ξ, η, ϕ, k, basis, even, oblate)
    end
    throw(ArgumentError("No hardcoded N reference implemented for $comp. Currently only :z is available."))
end

"""
    compare_N_generic_hardcoded(comp::Symbol, m, c, ξ, η, ϕ, k, basis;
        even=true, oblate=false, atol=1e-12)

Compare block-based `N_from_mode(comp,...)` against available hardcoded references.
Currently, hardcoded reference exists only for `comp == :z` (via `Nᶻ`).
"""
function compare_N_generic_hardcoded(comp::Symbol, m, c, ξ, η, ϕ, k, basis;
        even::Bool=true, oblate::Bool=false, atol=1e-12)
    gNξ, gNη, gNϕ = N_from_mode(comp, m, c, ξ, η, ϕ, k, basis; even=even, oblate=oblate)
    hNξ, hNη, hNϕ = _hardcoded_N_from_symbol(comp, m, c, ξ, η, ϕ, k, basis, even, oblate)

    Δξ = gNξ - hNξ
    Δη = gNη - hNη
    Δϕ = gNϕ - hNϕ

    relξ = abs(Δξ) / max(abs(gNξ), abs(hNξ), atol)
    relη = abs(Δη) / max(abs(gNη), abs(hNη), atol)
    relϕ = abs(Δϕ) / max(abs(gNϕ), abs(hNϕ), atol)

    return (
        generic   = (ξ=gNξ, η=gNη, ϕ=gNϕ),
        hardcoded = (ξ=hNξ, η=hNη, ϕ=hNϕ),
        abs_error = (ξ=abs(Δξ), η=abs(Δη), ϕ=abs(Δϕ)),
        rel_error = (ξ=relξ, η=relη, ϕ=relϕ),
    )
end

compare_Nz_generic_hardcoded(m, c, ξ, η, ϕ, k, basis; even::Bool=true, oblate::Bool=false, atol=1e-12) =
    compare_N_generic_hardcoded(:z, m, c, ξ, η, ϕ, k, basis; even=even, oblate=oblate, atol=atol)

"""
    compare_N_alpha_corr_sweep(comp::Symbol, m, c, η, ϕ, k, basis;
        even=true, oblate=false, ξs=nothing)

Sweep in ξ and report `(alpha, corr, err)` per spheroidal component.
"""
function compare_N_alpha_corr_sweep(comp::Symbol, m, c, η, ϕ, k, basis;
        even::Bool=true, oblate::Bool=false, ξs=nothing)

    ξs === nothing && (ξs = oblate ? range(0.05, 3.0, length=200) : range(1.05, 3.0, length=200))

    hξ = ComplexF64[]; hη = ComplexF64[]; hϕ = ComplexF64[]
    gξ = ComplexF64[]; gη = ComplexF64[]; gϕ = ComplexF64[]

    for ξ in ξs
        r = compare_N_generic_hardcoded(comp, m, c, ξ, η, ϕ, k, basis; even=even, oblate=oblate)
        push!(hξ, r.hardcoded.ξ); push!(hη, r.hardcoded.η); push!(hϕ, r.hardcoded.ϕ)
        push!(gξ, r.generic.ξ);   push!(gη, r.generic.η);   push!(gϕ, r.generic.ϕ)
    end

    return (
        ξ = (alpha = fit_alpha(hξ, gξ), corr = cosine_sim(hξ, gξ), err = shape_err(hξ, gξ)),
        η = (alpha = fit_alpha(hη, gη), corr = cosine_sim(hη, gη), err = shape_err(hη, gη)),
        ϕ = (alpha = fit_alpha(hϕ, gϕ), corr = cosine_sim(hϕ, gϕ), err = shape_err(hϕ, gϕ)),
    )
end

compare_Nz_alpha_corr_sweep(m, c, η, ϕ, k, basis; even::Bool=true, oblate::Bool=false, ξs=nothing) =
    compare_N_alpha_corr_sweep(:z, m, c, η, ϕ, k, basis; even=even, oblate=oblate, ξs=ξs)

