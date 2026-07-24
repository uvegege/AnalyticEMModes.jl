function check_confocal_elliptic_coax(focal_distance, inner_xi, outer_xi)
    focal_distance > 0 || throw(ArgumentError("focal_distance must be positive"))
    0 <= inner_xi < outer_xi || throw(ArgumentError("elliptic coax requires 0 <= inner_xi < outer_xi"))
    return nothing
end

"""
    EllipticCoax(inner_a, inner_b, outer_a, outer_b; atol=1e-10, rtol=1e-10)

Return `(focal_distance, inner_xi, outer_xi)` for a confocal elliptic coaxial
waveguide built from semi-axes. The two ellipses must be confocal, i.e.
`a^2 - b^2` must be equal for both boundaries.
"""
function EllipticCoax(inner_a, inner_b, outer_a, outer_b; atol=1.0e-10, rtol=1.0e-10)
    inner_a > inner_b > 0 || throw(ArgumentError("inner ellipse requires inner_a > inner_b > 0"))
    outer_a > outer_b > 0 || throw(ArgumentError("outer ellipse requires outer_a > outer_b > 0"))
    outer_a > inner_a && outer_b > inner_b || throw(ArgumentError("outer ellipse must contain the inner ellipse"))
    c2_inner = inner_a^2 - inner_b^2
    c2_outer = outer_a^2 - outer_b^2
    isapprox(c2_inner, c2_outer; atol=atol, rtol=rtol) || throw(ArgumentError("elliptic coax boundaries must be confocal: a^2 - b^2 must match"))
    c = sqrt((c2_inner + c2_outer) / 2)
    return c, acosh(inner_a / c), acosh(outer_a / c)
end

inner_axes_elliptic_coax(focal_distance, inner_xi) = (focal_distance * cosh(inner_xi), focal_distance * sinh(inner_xi))

outer_axes_elliptic_coax(focal_distance, outer_xi) = (focal_distance * cosh(outer_xi), focal_distance * sinh(outer_xi))

elliptic_coax_kind(kind) = kind == :TE || kind == :TM ? kind : throw(ArgumentError("elliptic coax mode kind must be :TE or :TM"))

function elliptic_coax_indices(m, n, even, kind)
    elliptic_coax_kind(kind)
    m >= 0 || throw(ArgumentError("Mathieu order m must be >= 0"))
    n >= 1 || throw(ArgumentError("radial index n must be >= 1"))
    (!even && m == 0) && throw(ArgumentError("odd Mathieu modes require m >= 1"))
    return nothing
end

"""
    ConfocalEllipticCoaxMode(focal_distance, inner_xi, outer_xi, m, n, even, kind)

Precomputed TE/TM mode for a confocal elliptic coaxial waveguide. This stores
the Mathieu parameter `q`, angular Fourier coefficients, radial coefficients,
and cutoff wavenumber, so field evaluation does not recompute Mathieu
characteristic values at each point.
"""
struct ConfocalEllipticCoaxMode{T, C}
    focal_distance::T
    inner_xi::T
    outer_xi::T
    m::Int
    n::Int
    even::Bool
    kind::Symbol
    q::T
    kc::T
    coeff::C
    radial_a::T
    radial_b::T
end

function elliptic_coax_coeff(m, q, even)
    a = even ? MathieuCharA(m, q) : MathieuCharB(m, q)
    return even ? mathieu_a_coeff(m, q, a, 100) : mathieu_b_coeff(m, q, a, 100)
end

function elliptic_coax_radial(kind_index, m, coeff, q, ξ, even)
    return even ? Mce_kernel(kind_index, m, coeff, q, ξ) : Mse_kernel(kind_index, m, coeff, q, ξ)
end

function elliptic_coax_angular(m, coeff, q, η, even)
    return even ? ce_kernel(m, coeff, q, η) : se_kernel(m, coeff, q, η)
end

function elliptic_coax_boundary_coefficients(m, q, coeff, ξ1, even, kind)
    R1, R1p = elliptic_coax_radial(1, m, coeff, q, ξ1, even)
    R2, R2p = elliptic_coax_radial(2, m, coeff, q, ξ1, even)
    return kind == :TE ? (1.0, -R1p / R2p) : (1.0, -R1 / R2)
end

function elliptic_coax_radial_value(mode::ConfocalEllipticCoaxMode, ξ)
    U1, U1p = elliptic_coax_radial(1, mode.m, mode.coeff, mode.q, ξ, mode.even)
    U2, U2p = elliptic_coax_radial(2, mode.m, mode.coeff, mode.q, ξ, mode.even)
    U = mode.radial_a * U1 + mode.radial_b * U2
    Up = mode.radial_a * U1p + mode.radial_b * U2p
    return U, Up
end

function elliptic_coax_determinant(q, focal_distance, inner_xi, outer_xi, m, even, kind)
    check_confocal_elliptic_coax(focal_distance, inner_xi, outer_xi)
    coeff = elliptic_coax_coeff(m, q, even)
    ξ1 = inner_xi
    ξ2 = outer_xi
    U11, U11p = elliptic_coax_radial(1, m, coeff, q, ξ1, even)
    U21, U21p = elliptic_coax_radial(2, m, coeff, q, ξ1, even)
    U12, U12p = elliptic_coax_radial(1, m, coeff, q, ξ2, even)
    U22, U22p = elliptic_coax_radial(2, m, coeff, q, ξ2, even)
    return kind == :TE ? U11p * U22p - U12p * U21p : U11 * U22 - U12 * U21
end

elliptic_coax_determinant_tuple(q, focal_distance, inner_xi, outer_xi, m, even, kind) = (elliptic_coax_determinant(q, focal_distance, inner_xi, outer_xi, m, even, kind), 0.0)

function q_elliptic_coax(focal_distance, inner_xi, outer_xi, m, n, even, kind)
    check_confocal_elliptic_coax(focal_distance, inner_xi, outer_xi)
    elliptic_coax_indices(m, n, even, kind)
    f(q) = elliptic_coax_determinant_tuple(q, focal_distance, inner_xi, outer_xi, m, even, kind)
    return find_n_zero_ewg(f, m, n, kind; inc_iters=900, max_iters=90, tol=1.0e-10)
end

"""
    kc_elliptic_coax(focal_distance, inner_xi, outer_xi, m, n, even, kind)

Cutoff wavenumber for TE/TM modes in a confocal elliptic coaxial waveguide.
For the TEM mode use `kc_elliptic_coax(focal_distance, inner_xi, outer_xi, :TEM)`.
"""
function kc_elliptic_coax(focal_distance, inner_xi, outer_xi, m, n, even, kind)
    q = q_elliptic_coax(focal_distance, inner_xi, outer_xi, m, n, even, kind)
    return 2 * sqrt(q) / focal_distance
end

kc_elliptic_coax(focal_distance, inner_xi, outer_xi, ::Val{:TEM}) = (check_confocal_elliptic_coax(focal_distance, inner_xi, outer_xi); 0.0)

kc_elliptic_coax(focal_distance, inner_xi, outer_xi, kind::Symbol) = kind == :TEM ? kc_elliptic_coax(focal_distance, inner_xi, outer_xi, Val(:TEM)) : throw(ArgumentError("only :TEM is supported by this method"))

function ConfocalEllipticCoaxMode(focal_distance, inner_xi, outer_xi, m, n, even, kind)
    check_confocal_elliptic_coax(focal_distance, inner_xi, outer_xi)
    elliptic_coax_indices(m, n, even, kind)
    q = q_elliptic_coax(focal_distance, inner_xi, outer_xi, m, n, even, kind)
    kc = 2 * sqrt(q) / focal_distance
    coeff = elliptic_coax_coeff(m, q, even)
    radial_a, radial_b = elliptic_coax_boundary_coefficients(m, q, coeff, inner_xi, even, kind)
    return ConfocalEllipticCoaxMode(focal_distance, inner_xi, outer_xi, m, n, even, kind, q, kc, coeff, radial_a, radial_b)
end

mode_label_elliptic_coax(mode::ConfocalEllipticCoaxMode) = mode.kind == :TE ? (mode.even ? :TEe : :TEo) : (mode.even ? :TMe : :TMo)

function elliptic_coax_modal_f(mode::ConfocalEllipticCoaxMode, ξ, η)
    U, Uξ = elliptic_coax_radial_value(mode, ξ)
    V, Vη = elliptic_coax_angular(mode.m, mode.coeff, mode.q, η, mode.even)
    return (Uξ * V, U * Vη, U * V)
end

ψ(mode::ConfocalEllipticCoaxMode, ξ, η) = elliptic_coax_modal_f(mode, ξ, η)[3]

∇ψ(mode::ConfocalEllipticCoaxMode, ξ, η) = elliptic_coax_modal_f(mode, ξ, η)[1:2]

function te_elliptic_coax_fields(mode::ConfocalEllipticCoaxMode, ξ, η, c_e, c_h)
    mode.kind == :TE || throw(ArgumentError("mode must be TE"))
    ∂ψᵢ, ∂ψⱼ, ψₖ = elliptic_coax_modal_f(mode, ξ, η)
    h = mode.focal_distance * sqrt(sinh(ξ)^2 + sin(η)^2)
    Eξ = -c_e / h * ∂ψⱼ
    Eη = +c_e / h * ∂ψᵢ
    Ez = zero(Eξ)
    Hξ = -c_h / h * ∂ψᵢ
    Hη = -c_h / h * ∂ψⱼ
    Hz = -im * ψₖ
    return (Eξ, Eη, Ez, Hξ, Hη, Hz)
end

function tm_elliptic_coax_fields(mode::ConfocalEllipticCoaxMode, ξ, η, c_e, c_h)
    mode.kind == :TM || throw(ArgumentError("mode must be TM"))
    ∂ψᵢ, ∂ψⱼ, ψₖ = elliptic_coax_modal_f(mode, ξ, η)
    h = mode.focal_distance * sqrt(sinh(ξ)^2 + sin(η)^2)
    Eξ = -c_e / h * ∂ψᵢ
    Eη = -c_e / h * ∂ψⱼ
    Ez = -im * ψₖ
    Hξ = +c_h / h * ∂ψⱼ
    Hη = -c_h / h * ∂ψᵢ
    Hz = zero(Hξ)
    return (Eξ, Eη, Ez, Hξ, Hη, Hz)
end

function te_elliptic_coax_fields(mode::ConfocalEllipticCoaxMode, ξ, η, f, μᵣ, εᵣ)
    β = phase_constant(mode.kc, f, μᵣ, εᵣ)
    c_e, c_h = te_coefficients(mode.kc, β, f, μᵣ, εᵣ)
    return te_elliptic_coax_fields(mode, ξ, η, c_e, c_h)
end

function tm_elliptic_coax_fields(mode::ConfocalEllipticCoaxMode, ξ, η, f, μᵣ, εᵣ)
    β = phase_constant(mode.kc, f, μᵣ, εᵣ)
    c_e, c_h = tm_coefficients(mode.kc, β, f, μᵣ, εᵣ)
    return tm_elliptic_coax_fields(mode, ξ, η, c_e, c_h)
end

te_elliptic_coax_fields(mode::ConfocalEllipticCoaxMode, ξ::AbstractArray, η::AbstractArray, f, μᵣ, εᵣ) = te_elliptic_coax_fields.(Ref(mode), ξ, η, f, μᵣ, εᵣ)

tm_elliptic_coax_fields(mode::ConfocalEllipticCoaxMode, ξ::AbstractArray, η::AbstractArray, f, μᵣ, εᵣ) = tm_elliptic_coax_fields.(Ref(mode), ξ, η, f, μᵣ, εᵣ)

function te_elliptic_coax_fields(ξ, η, focal_distance, inner_xi, outer_xi, m, n, even, f, μᵣ, εᵣ)
    return te_elliptic_coax_fields(ConfocalEllipticCoaxMode(focal_distance, inner_xi, outer_xi, m, n, even, :TE), ξ, η, f, μᵣ, εᵣ)
end

function tm_elliptic_coax_fields(ξ, η, focal_distance, inner_xi, outer_xi, m, n, even, f, μᵣ, εᵣ)
    return tm_elliptic_coax_fields(ConfocalEllipticCoaxMode(focal_distance, inner_xi, outer_xi, m, n, even, :TM), ξ, η, f, μᵣ, εᵣ)
end

function tem_elliptic_coax_fields(ξ, η, focal_distance, inner_xi, outer_xi, μᵣ, εᵣ)
    check_confocal_elliptic_coax(focal_distance, inner_xi, outer_xi)
    ε = εᵣ * _εₒ
    μ = μᵣ * _μₒ
    η₀ = sqrt(μ / ε)
    Δξ = outer_xi - inner_xi
    Z0 = η₀ * Δξ / (2 * π)
    h = focal_distance * sqrt(sinh(ξ)^2 + sin(η)^2)
    Eξ = sqrt(2 * Z0) / (h * Δξ)
    Eη = zero(Eξ)
    Ez = zero(Eξ)
    Hξ = zero(Eξ)
    Hη = Eξ / η₀
    Hz = zero(Eξ)
    return (Eξ, Eη, Ez, Hξ, Hη, Hz)
end

function elliptic_coax_modal_power_1d(mode::ConfocalEllipticCoaxMode, β, f, μᵣ, εᵣ; Nη=4096, rtol=1.0e-10, atol=0.0, maxdepth=25)
    A_v, A_vp = Av_Avp_from_kernels(mode.m, mode.q, mode.coeff; even=mode.even, Nη=Nη)
    fU(ξ) = abs2(elliptic_coax_radial_value(mode, ξ)[1])
    fUp(ξ) = abs2(elliptic_coax_radial_value(mode, ξ)[2])
    Iu = quad_asimpson(fU, mode.inner_xi, mode.outer_xi; rtol=rtol, atol=atol, maxdepth=maxdepth)
    Iup = quad_asimpson(fUp, mode.inner_xi, mode.outer_xi; rtol=rtol, atol=atol, maxdepth=maxdepth)
    ω = 2 * π * f
    factor = mode.kind == :TE ? ω * μᵣ * _μₒ * β / mode.kc^4 : ω * εᵣ * _εₒ * β / mode.kc^4
    return 0.5 * factor * (A_v * Iup + A_vp * Iu)
end

te_normalization_elliptic_coax(mode::ConfocalEllipticCoaxMode, β, f, μᵣ, εᵣ; kwargs...) = mode.kind == :TE ? sqrt(1 / elliptic_coax_modal_power_1d(mode, β, f, μᵣ, εᵣ; kwargs...)) : throw(ArgumentError("mode must be TE"))

tm_normalization_elliptic_coax(mode::ConfocalEllipticCoaxMode, β, f, μᵣ, εᵣ; kwargs...) = mode.kind == :TM ? sqrt(1 / elliptic_coax_modal_power_1d(mode, β, f, μᵣ, εᵣ; kwargs...)) : throw(ArgumentError("mode must be TM"))

function te_normalization_elliptic_coax(focal_distance, inner_xi, outer_xi, m, n, even, β, f, μᵣ, εᵣ; kwargs...)
    mode = ConfocalEllipticCoaxMode(focal_distance, inner_xi, outer_xi, m, n, even, :TE)
    return te_normalization_elliptic_coax(mode, β, f, μᵣ, εᵣ; kwargs...)
end

function tm_normalization_elliptic_coax(focal_distance, inner_xi, outer_xi, m, n, even, β, f, μᵣ, εᵣ; kwargs...)
    mode = ConfocalEllipticCoaxMode(focal_distance, inner_xi, outer_xi, m, n, even, :TM)
    return tm_normalization_elliptic_coax(mode, β, f, μᵣ, εᵣ; kwargs...)
end

function first_n_modes_elliptic_coax(N, focal_distance, inner_xi, outer_xi)
    check_confocal_elliptic_coax(focal_distance, inner_xi, outer_xi)
    N <= 0 && return Tuple{Symbol, Int, Int, Float64}[]
    modes = Tuple{Symbol, Int, Int, Float64}[]
    mmax = max(2, ceil(Int, sqrt(N)) + 2)
    nmax = max(2, ceil(Int, sqrt(N)) + 2)
    while true
        empty!(modes)
        push!(modes, (:TEM, 0, 0, 0.0))
        for m in 0:mmax, n in 1:nmax
            for kind in (:TE, :TM)
                mode = ConfocalEllipticCoaxMode(focal_distance, inner_xi, outer_xi, m, n, true, kind)
                push!(modes, (mode_label_elliptic_coax(mode), m, n, mode.kc))
            end
            if m > 0
                for kind in (:TE, :TM)
                    mode = ConfocalEllipticCoaxMode(focal_distance, inner_xi, outer_xi, m, n, false, kind)
                    push!(modes, (mode_label_elliptic_coax(mode), m, n, mode.kc))
                end
            end
        end
        sort!(modes, by=x -> (x[4], x[1], x[2], x[3]))
        if length(modes) >= N
            selected = modes[1:N]
            touches_mmax = any(mode -> mode[2] == mmax, selected)
            touches_nmax = any(mode -> mode[3] == nmax, selected)
            !(touches_mmax || touches_nmax) && break
            touches_mmax && (mmax *= 2)
            touches_nmax && (nmax *= 2)
        else
            mmax *= 2
            nmax *= 2
        end
    end
    return modes[1:N]
end
