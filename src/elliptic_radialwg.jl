function check_elliptic_radial(focal_distance, height)
    focal_distance > 0 || throw(ArgumentError("focal_distance must be positive"))
    height > 0 || throw(ArgumentError("height must be positive"))
    return nothing
end

phase_constant_elliptic_radial(height, n) = n * π / height

cutoff_frequency_elliptic_radial(height, n, μᵣ, εᵣ) = cutoff_frequency(phase_constant_elliptic_radial(height, n), μᵣ, εᵣ)

"""
    kt_elliptic_radial(focal_distance, height, n, f, μᵣ, εᵣ)

Transverse elliptic-radial propagation constant. This is the analogue of
`kc_radial`: `k_t = sqrt(k^2 - (n*pi/height)^2)`.
"""
function kt_elliptic_radial(focal_distance, height, n, f, μᵣ, εᵣ)
    check_elliptic_radial(focal_distance, height)
    βz = phase_constant_elliptic_radial(height, n)
    μ = μᵣ * _μₒ
    ε = εᵣ * _εₒ
    k = 2 * π * f * sqrt(μ * ε)
    value = k^2 - βz^2
    value > 0 || throw(ArgumentError("elliptic radial open modes require f above cutoff"))
    return sqrt(value)
end

function elliptic_radial_indices(m, n, even, kind)
    kind == :TE || kind == :TM || throw(ArgumentError("elliptic radial mode kind must be :TE or :TM"))
    m >= 0 || throw(ArgumentError("Mathieu order m must be >= 0"))
    n >= 0 || throw(ArgumentError("plate index n must be >= 0"))
    kind == :TE && n == 0 && throw(ArgumentError("TE elliptic radial modes require n >= 1"))
    !even && m == 0 && throw(ArgumentError("odd Mathieu modes require m >= 1"))
    return nothing
end

"""
    EllipticRadialMode(focal_distance, height, m, n, even, kind, f, μᵣ, εᵣ; outgoing=true)

Precomputed open elliptic radial mode. It stores `k_t`, `q`, angular Mathieu
coefficients, and the outgoing/incoming radial Mathieu-Hankel sign, so field
evaluation at many points does not recompute Mathieu characteristic data.
"""
struct EllipticRadialMode{T, Q, C}
    focal_distance::T
    height::T
    m::Int
    n::Int
    even::Bool
    kind::Symbol
    kt::Q
    q::Q
    coeff::C
    hankel_sign::Int
end

function elliptic_radial_coeff(m, q, even)
    a = even ? MathieuCharA(m, q) : MathieuCharB(m, q)
    return even ? mathieu_a_coeff(m, q, a, 100) : mathieu_b_coeff(m, q, a, 100)
end

function EllipticRadialMode(focal_distance, height, m, n, even, kind, f, μᵣ, εᵣ; outgoing=true)
    check_elliptic_radial(focal_distance, height)
    elliptic_radial_indices(m, n, even, kind)
    kt = kt_elliptic_radial(focal_distance, height, n, f, μᵣ, εᵣ)
    q = (kt * focal_distance)^2 / 4
    abs(q) > 0 || throw(ArgumentError("elliptic radial Mathieu-Hankel fields are singular at exact cutoff"))
    coeff = elliptic_radial_coeff(m, q, even)
    hankel_sign = outgoing ? 1 : -1
    return EllipticRadialMode(focal_distance, height, m, n, even, kind, kt, q, coeff, hankel_sign)
end

function elliptic_radial_angular(mode::EllipticRadialMode, η)
    return mode.even ? ce_kernel(mode.m, mode.coeff, mode.q, η) : se_kernel(mode.m, mode.coeff, mode.q, η)
end

function elliptic_radial_radial(mode::EllipticRadialMode, ξ)
    U1, U1ξ = mode.even ? Mce_kernel(1, mode.m, mode.coeff, mode.q, ξ) : Mse_kernel(1, mode.m, mode.coeff, mode.q, ξ)
    U2, U2ξ = mode.even ? Mce_kernel(2, mode.m, mode.coeff, mode.q, ξ) : Mse_kernel(2, mode.m, mode.coeff, mode.q, ξ)
    s = im * mode.hankel_sign
    return U1 + s * U2, U1ξ + s * U2ξ
end

function elliptic_radial_modal_f(mode::EllipticRadialMode, ξ, η)
    U, Uξ = elliptic_radial_radial(mode, ξ)
    V, Vη = elliptic_radial_angular(mode, η)
    return (Uξ * V, U * Vη, U * V)
end

ψ(mode::EllipticRadialMode, ξ, η) = elliptic_radial_modal_f(mode, ξ, η)[3]

∇ψ(mode::EllipticRadialMode, ξ, η) = elliptic_radial_modal_f(mode, ξ, η)[1:2]

function te_elliptic_radial_fields(mode::EllipticRadialMode, ξ, η, z, c_e, c_h)
    mode.kind == :TE || throw(ArgumentError("mode must be TE"))
    ∂ψᵢ, ∂ψⱼ, ψₖ = elliptic_radial_modal_f(mode, ξ, η)
    βz = phase_constant_elliptic_radial(mode.height, mode.n)
    snz, cnz = sincos(βz * z)
    h = mode.focal_distance * sqrt(sinh(ξ)^2 + sin(η)^2)
    Eξ = -c_e / h * ∂ψⱼ * snz
    Eη = +c_e / h * ∂ψᵢ * snz
    Ez = zero(Eξ)
    Hξ = -im * c_h / h * ∂ψᵢ * cnz
    Hη = -im * c_h / h * ∂ψⱼ * cnz
    Hz = -im * ψₖ * snz
    return (Eξ, Eη, Ez, Hξ, Hη, Hz)
end

function tm_elliptic_radial_fields(mode::EllipticRadialMode, ξ, η, z, c_e, c_h)
    mode.kind == :TM || throw(ArgumentError("mode must be TM"))
    ∂ψᵢ, ∂ψⱼ, ψₖ = elliptic_radial_modal_f(mode, ξ, η)
    βz = phase_constant_elliptic_radial(mode.height, mode.n)
    snz, cnz = sincos(βz * z)
    h = mode.focal_distance * sqrt(sinh(ξ)^2 + sin(η)^2)
    Eξ = +im * c_e / h * ∂ψᵢ * snz
    Eη = +im * c_e / h * ∂ψⱼ * snz
    Ez = -im * ψₖ * cnz
    Hξ = +c_h / h * ∂ψⱼ * cnz
    Hη = -c_h / h * ∂ψᵢ * cnz
    Hz = zero(Eξ)
    return (Eξ, Eη, Ez, Hξ, Hη, Hz)
end

function te_elliptic_radial_fields(mode::EllipticRadialMode, ξ, η, z, f, μᵣ, εᵣ)
    βz = phase_constant_elliptic_radial(mode.height, mode.n)
    c_e, c_h = te_coefficients(mode.kt, βz, f, μᵣ, εᵣ)
    return te_elliptic_radial_fields(mode, ξ, η, z, c_e, c_h)
end

function tm_elliptic_radial_fields(mode::EllipticRadialMode, ξ, η, z, f, μᵣ, εᵣ)
    βz = phase_constant_elliptic_radial(mode.height, mode.n)
    c_e, c_h = tm_coefficients(mode.kt, βz, f, μᵣ, εᵣ)
    return tm_elliptic_radial_fields(mode, ξ, η, z, c_e, c_h)
end

te_elliptic_radial_fields(mode::EllipticRadialMode, ξ::AbstractArray, η::AbstractArray, z::AbstractArray, f, μᵣ, εᵣ) = te_elliptic_radial_fields.(Ref(mode), ξ, η, z, f, μᵣ, εᵣ)

tm_elliptic_radial_fields(mode::EllipticRadialMode, ξ::AbstractArray, η::AbstractArray, z::AbstractArray, f, μᵣ, εᵣ) = tm_elliptic_radial_fields.(Ref(mode), ξ, η, z, f, μᵣ, εᵣ)

function te_elliptic_radial_fields(ξ, η, z, focal_distance, height, m, n, even, f, μᵣ, εᵣ; outgoing=true)
    mode = EllipticRadialMode(focal_distance, height, m, n, even, :TE, f, μᵣ, εᵣ; outgoing=outgoing)
    return te_elliptic_radial_fields(mode, ξ, η, z, f, μᵣ, εᵣ)
end

function tm_elliptic_radial_fields(ξ, η, z, focal_distance, height, m, n, even, f, μᵣ, εᵣ; outgoing=true)
    mode = EllipticRadialMode(focal_distance, height, m, n, even, :TM, f, μᵣ, εᵣ; outgoing=outgoing)
    return tm_elliptic_radial_fields(mode, ξ, η, z, f, μᵣ, εᵣ)
end

"""
    elliptic_radial_modal_power_1d(mode, ξ, f, μᵣ, εᵣ; Nη=4096)

Time-averaged power crossing the open elliptic-radial surface `ξ = const`
for a mode with unit longitudinal amplitude.

The integral is the elliptic-coordinate analogue of the radial waveguide
normalization at `r = const`. Since the flux is through a constant-`ξ`
surface, the metric factor cancels between the transverse field and the
surface element, and the result only needs the angular integral
`A_v = ∫₀²π |V(η)|² dη` and the radial Mathieu-Hankel Wronskian term.

# Arguments
- `mode`: Precomputed elliptic-radial mode
- `ξ`: Elliptic radial coordinate of the normalization surface
- `f`: Frequency in Hz
- `μᵣ`: Relative permeability
- `εᵣ`: Relative permittivity

# Keyword Arguments
- `Nη`: Number of angular quadrature points for `A_v` (default 4096)
"""
function elliptic_radial_modal_power_1d(mode::EllipticRadialMode, ξ, f, μᵣ, εᵣ; Nη::Int=4096)
    ω = 2π * f
    σ = mode.kind == :TE ? μᵣ * _μₒ :
        mode.kind == :TM ? εᵣ * _εₒ :
        throw(ArgumentError("mode must be TE or TM"))
    U, Uξ = elliptic_radial_radial(mode, ξ)
    A_v, _ = Av_Avp_from_kernels(mode.m, mode.q, mode.coeff; even=mode.even, Nη=Nη)
    P = 1/2 * ω * σ / mode.kt^2 * mode.height / 2 * A_v * real(im * conj(U) * Uξ)
    return abs(P)
end

"""
    te_normalization_elliptic_radial(mode, ξ, f, μᵣ, εᵣ; kwargs...)

Normalization factor for TE modes in an elliptic radial guide to achieve unit
power through the open surface `ξ = const`.

Returns `√(1/P)`, where `P` is the time-averaged power computed by
[`elliptic_radial_modal_power_1d`].
"""
function te_normalization_elliptic_radial(mode::EllipticRadialMode, ξ, f, μᵣ, εᵣ; kwargs...)
    mode.kind == :TE || throw(ArgumentError("mode must be TE"))
    return sqrt(1 / elliptic_radial_modal_power_1d(mode, ξ, f, μᵣ, εᵣ; kwargs...))
end

"""
    tm_normalization_elliptic_radial(mode, ξ, f, μᵣ, εᵣ; kwargs...)

Normalization factor for TM modes in an elliptic radial guide to achieve unit
power through the open surface `ξ = const`.

Returns `√(1/P)`, where `P` is the time-averaged power computed by
[`elliptic_radial_modal_power_1d`].
"""
function tm_normalization_elliptic_radial(mode::EllipticRadialMode, ξ, f, μᵣ, εᵣ; kwargs...)
    mode.kind == :TM || throw(ArgumentError("mode must be TM"))
    return sqrt(1 / elliptic_radial_modal_power_1d(mode, ξ, f, μᵣ, εᵣ; kwargs...))
end

function te_normalization_elliptic_radial(ξ, focal_distance, height, m, n, even, f, μᵣ, εᵣ; outgoing=true, kwargs...)
    mode = EllipticRadialMode(focal_distance, height, m, n, even, :TE, f, μᵣ, εᵣ; outgoing=outgoing)
    return te_normalization_elliptic_radial(mode, ξ, f, μᵣ, εᵣ; kwargs...)
end

function tm_normalization_elliptic_radial(ξ, focal_distance, height, m, n, even, f, μᵣ, εᵣ; outgoing=true, kwargs...)
    mode = EllipticRadialMode(focal_distance, height, m, n, even, :TM, f, μᵣ, εᵣ; outgoing=outgoing)
    return tm_normalization_elliptic_radial(mode, ξ, f, μᵣ, εᵣ; kwargs...)
end

"""
    first_n_modes_elliptic_radial(N, height, μᵣ=1.0, εᵣ=1.0)

Return the first `N` cutoff entries of the open elliptic radial guide. Since
the cutoff depends only on the plate index `n`, each entry uses `m = 0` as a
representative of the infinitely many Mathieu transverse shapes sharing that
same cutoff.
"""
function first_n_modes_elliptic_radial(N, height, μᵣ=1.0, εᵣ=1.0)
    height > 0 || throw(ArgumentError("height must be positive"))
    N <= 0 && return Tuple{Symbol, Int, Int, Float64}[]
    modes = Tuple{Symbol, Int, Int, Float64}[]
    n = 0
    while length(modes) < N
        n > 0 && push!(modes, (:TE, 0, n, cutoff_frequency_elliptic_radial(height, n, μᵣ, εᵣ)))
        push!(modes, (:TM, 0, n, cutoff_frequency_elliptic_radial(height, n, μᵣ, εᵣ)))
        n += 1
    end
    return modes[1:N]
end
