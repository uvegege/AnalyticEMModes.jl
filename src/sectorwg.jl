sector_order(p, ϕ0) = p * π / ϕ0

function check_sector_kind(kind)
    if !(kind == :TE || kind == :TM)
        throw(ArgumentError("sector mode kind must be :TE or :TM"))
    end
    return nothing
end

function check_sector_indices(p, n, kind)
    check_sector_kind(kind)
    n >= 1 || throw(ArgumentError("radial index n must be >= 1"))
    p >= 0 || throw(ArgumentError("angular index p must be >= 0"))
    kind == :TM && p == 0 && throw(ArgumentError("sector TM modes require p >= 1"))
    return nothing
end

function check_sector_angle(ϕ0)
    0 < ϕ0 < 2 * π || throw(ArgumentError("sector angle must satisfy 0 < angle < 2*pi"))
    return nothing
end

function bisect_root(f, a, b; tol=1.0e-12, max_iters=100)
    fa = f(a)
    fb = f(b)
    iszero(fa) && return a
    iszero(fb) && return b
    sign(fa) == sign(fb) && throw(ArgumentError("root is not bracketed"))
    left = a
    right = b
    for _ in 1:max_iters
        mid = (left + right) / 2
        fm = f(mid)
        abs(fm) < tol && return mid
        if sign(fm) == sign(fa)
            left = mid
            fa = fm
        else
            right = mid
        end
        abs(right - left) < tol * max(1, abs(mid)) && return (left + right) / 2
    end
    return (left + right) / 2
end

function nth_positive_root(f, n; x0=1.0e-7, step=π / 32, tol=1.0e-12)
    roots = Float64[]
    x_left = x0
    f_left = f(x_left)
    x = x_left + step
    while length(roots) < n
        f_right = f(x)
        if isfinite(f_left) && isfinite(f_right) && !iszero(sign(f_left)) && !iszero(sign(f_right)) && sign(f_left) != sign(f_right)
            root = bisect_root(f, x_left, x; tol=tol)
            if isempty(roots) || abs(root - roots[end]) > sqrt(tol) * max(1, root)
                push!(roots, root)
            end
        end
        x_left = x
        f_left = f_right
        x += step
        x > x0 + (n + 40) * π + 4 * step && length(roots) < n && (step *= 0.5)
        x > x0 + (n + 200) * π && error("could not locate root $n")
    end
    return roots[n]
end

root_search_start(ν) = ν > 4 ? ν - 2 : (ν == 0 ? 1.0e-5 : 1.0e-6)

function besselj_zero_real_order(ν, n)
    return nth_positive_root(x -> besselj(ν, x), n; x0=root_search_start(ν))
end

function besselj_prime_zero_real_order(ν, n)
    x0 = ν == 0 ? root_search_start(ν) : max(0.01, root_search_start(ν))
    return nth_positive_root(x -> besselj_prime(ν, x), n; x0=x0)
end

function mcmahon_besselj_zero_seed(ν, n)
    a = π * (n + ν / 2 - 1 / 4)
    μ = 4 * ν^2
    return a - (μ - 1) / (8 * a) - (4 * μ - 1) * (7 * μ - 31) / (3 * (8 * a)^3) -
           32 * (μ - 1) * (83 * μ^2 - 982 * μ + 3779) / (15 * (8 * a)^5) -
           64 * (μ - 1) * (6949 * μ^3 - 153855 * μ^2 + 1585743 * μ - 6277237) / (105 * (8 * a)^7)
end

function mcmahon_besselj_prime_zero_seed(ν, n)
    ν == 0 && return mcmahon_besselj_zero_seed(1.0, n)
    b = π * (n + ν / 2 - 3 / 4)
    μ = 4 * ν^2
    return b - (μ + 3) / (8 * b) - 4 * (7 * μ^2 + 82 * μ - 9) / (3 * (8 * b)^3) -
           32 * (83 * μ^3 + 2075 * μ^2 - 3039 * μ + 3537) / (15 * (8 * b)^5) -
           64 * (6949 * μ^4 + 296492 * μ^3 - 1248002 * μ^2 + 7414380 * μ - 5853627) / (105 * (8 * b)^7)
end

const AIRY_ZEROS_SECTOR = (-2.338107410459767, -4.087949444130971, -5.520559828095551, -6.786708090071759, -7.944133587120853, -9.022650853340980, -10.04017434155809, -11.00852430373326)

const AIRY_PRIME_ZEROS_SECTOR = (-1.018792971647471, -3.248197582179837, -4.820099211178736, -6.163307355639487, -7.372177255047770, -8.488486734019722, -9.535449052433547, -10.52766039695740)

function airy_zero_sector(n)
    n <= length(AIRY_ZEROS_SECTOR) && return AIRY_ZEROS_SECTOR[n]
    t = 3 * π * (4 * n - 1) / 8
    return -t^(2 / 3) * (1 + 5 / (48 * t^2))
end

function airy_prime_zero_sector(n)
    n <= length(AIRY_PRIME_ZEROS_SECTOR) && return AIRY_PRIME_ZEROS_SECTOR[n]
    u = 3 * π * (4 * n - 3) / 8
    return -u^(2 / 3) * (1 - 7 / (48 * u^2))
end

function olver_delta_sector(kind, n)
    return 2.0^(-1 / 3) * -(kind == :TE ? airy_prime_zero_sector(n) : airy_zero_sector(n))
end

function olver_bessel_zero_seed(kind, ν, n)
    δ = olver_delta_sector(kind, n)
    c = kind == :TM ? 3 * δ^2 / 10 : 3 * δ^2 / 10 - 1 / (10 * δ)
    return ν + δ * ν^(1 / 3) + c * ν^(-1 / 3)
end

use_olver_sector_seed(ν, n) = ν >= 3 * n + 4

function sector_newton_clamp(kind, ν, n, use_olver)
    if use_olver
        gap = (olver_delta_sector(kind, n + 1) - olver_delta_sector(kind, n)) * ν^(1 / 3)
        return max(π, gap) / 2
    end
    return π / 2
end

function newton_bessel_zero(kind, ν, x0; clamp=π / 2, tol=1.0e-13, max_iters=40)
    x = x0
    for _ in 1:max_iters
        if kind == :TE
            f = besselj_prime(ν, x)
            df = besselj_doubleprime(ν, x)
        else
            f = besselj(ν, x)
            df = besselj_prime(ν, x)
        end
        Δ = f / df
        abs(Δ) > clamp && (Δ = sign(Δ) * clamp)
        x -= Δ
        abs(Δ) <= tol * max(1, abs(x)) && return x
    end
    return x
end

function sector_bessel_zero_residual(kind, ν, x)
    return kind == :TE ? abs(besselj_prime(ν, x)) : abs(besselj(ν, x))
end

function fast_bessel_zero_real_order(kind, ν, n)
    kind == :TE && ν == 0 && return fast_bessel_zero_real_order(:TM, 1.0, n)
    use_olver = use_olver_sector_seed(ν, n)
    x0 = use_olver ? olver_bessel_zero_seed(kind, ν, n) : (kind == :TE ? mcmahon_besselj_prime_zero_seed(ν, n) : mcmahon_besselj_zero_seed(ν, n))
    x = newton_bessel_zero(kind, ν, x0; clamp=sector_newton_clamp(kind, ν, n, use_olver))
    if isfinite(x) && x > 0 && sector_bessel_zero_residual(kind, ν, x) < 1.0e-10
        return x
    end
    return kind == :TE ? besselj_prime_zero_real_order(ν, n) : besselj_zero_real_order(ν, n)
end

"""
    kc_sector(radius, angle, p, n, kind)

Cutoff wavenumber of a circular-sector waveguide mode. `kind` is `:TE` or `:TM`.
"""
function kc_sector(radius, angle, p, n, kind)
    check_sector_angle(angle)
    check_sector_indices(p, n, kind)
    ν = sector_order(p, angle)
    x = fast_bessel_zero_real_order(kind, ν, n)
    return x / radius
end

characteristic_annular_sector_tm(ν, kc, b, a) = characteristic_coax_equation_tm(ν, kc, a, b)
characteristic_annular_sector_te(ν, kc, b, a) = characteristic_coax_equation_te(ν, kc, a, b)

function annular_sector_zero(ν, b, a, n, kind)
    f = kind == :TE ? x -> characteristic_annular_sector_te(ν, x, b, a) : x -> characteristic_annular_sector_tm(ν, x, b, a)
    return nth_positive_root(f, n; x0=1.0e-7 / a, step=π / (32 * (a - b)))
end

"""
    kc_annular_sector(inner_radius, outer_radius, angle, p, n, kind)

Cutoff wavenumber of an annular-sector waveguide mode. `kind` is `:TE` or `:TM`.
"""
function kc_annular_sector(inner_radius, outer_radius, angle, p, n, kind)
    check_sector_angle(angle)
    0 < inner_radius < outer_radius || throw(ArgumentError("annular sector requires 0 < inner_radius < outer_radius"))
    check_sector_indices(p, n, kind)
    ν = sector_order(p, angle)
    return annular_sector_zero(ν, inner_radius, outer_radius, n, kind)
end

function sector_angular(ϕ, ν, kind)
    if kind == :TE
        s, c = sincos(ν * ϕ)
        return c, -ν * s
    else
        s, c = sincos(ν * ϕ)
        return s, ν * c
    end
end

function nu_over_r_besselj(ν, kc, r)
    if iszero(r)
        iszero(ν) && return 0.0
        isapprox(ν, 1.0; atol=1.0e-14, rtol=1.0e-14) && return kc / 2
        ν > 1 && return 0.0
        return Inf
    end
    return ν * besselj(ν, kc * r) / r
end

function sector_modal_f(r, ϕ, radius, angle, p, n, kind)
    kc = kc_sector(radius, angle, p, n, kind)
    ν = sector_order(p, angle)
    Φ, Φp = sector_angular(ϕ, ν, kind)
    R = besselj(ν, kc * r)
    Rp = kc * besselj_prime(ν, kc * r)
    ∂ψᵣ = Rp * Φ
    ∂ψᵩ = iszero(Φp) ? zero(R) : R * Φp / r
    if iszero(r) && !iszero(Φp)
        ∂ψᵩ = nu_over_r_besselj(ν, kc, r) * (Φp / ν)
    end
    ψₖ = R * Φ
    return (∂ψᵣ, ∂ψᵩ, ψₖ)
end

annular_sector_boundary_coeff_tm(ν, b, kc) = coax_boundary_coeff_tm(ν, b, kc)

annular_sector_boundary_coeff_te(ν, b, kc) = coax_boundary_coeff_te(ν, b, kc)

function annular_sector_modal_f(r, ϕ, inner_radius, outer_radius, angle, p, n, kind)
    kc = kc_annular_sector(inner_radius, outer_radius, angle, p, n, kind)
    ν = sector_order(p, angle)
    Am, Bm = kind == :TE ? annular_sector_boundary_coeff_te(ν, inner_radius, kc) : annular_sector_boundary_coeff_tm(ν, inner_radius, kc)
    Φ, Φp = sector_angular(ϕ, ν, kind)
    R = Am * besselj(ν, kc * r) + Bm * bessely(ν, kc * r)
    Rp = kc * (Am * besselj_prime(ν, kc * r) + Bm * bessely_prime(ν, kc * r))
    return (Rp * Φ, R * Φp / r, R * Φ)
end

function sector_fields_from_modal(kind, ∂ψᵣ, ∂ψᵩ, ψₖ, c_e, c_h)
    if kind == :TE
        Er = -c_e * ∂ψᵩ
        Eϕ = +c_e * ∂ψᵣ
        Ez = zero(Er)
        Hr = -c_h * ∂ψᵣ
        Hϕ = -c_h * ∂ψᵩ
        Hz = -im * ψₖ
        return (Er, Eϕ, Ez, Hr, Hϕ, Hz)
    else
        Er = -c_e * ∂ψᵣ
        Eϕ = -c_e * ∂ψᵩ
        Ez = -im * ψₖ
        Hr = +c_h * ∂ψᵩ
        Hϕ = -c_h * ∂ψᵣ
        Hz = zero(Ez)
        return (Er, Eϕ, Ez, Hr, Hϕ, Hz)
    end
end

function sector_fields(kind, r, ϕ, kc, modal, f, μᵣ, εᵣ)
    β = phase_constant(kc, f, μᵣ, εᵣ)
    c_e, c_h = kind == :TE ? te_coefficients(kc, β, f, μᵣ, εᵣ) : tm_coefficients(kc, β, f, μᵣ, εᵣ)
    return sector_fields_from_modal(kind, modal..., c_e, c_h)
end

te_sector_fields(r, ϕ, radius, angle, p, n, f, μᵣ, εᵣ) = sector_fields(:TE, r, ϕ, kc_sector(radius, angle, p, n, :TE), sector_modal_f(r, ϕ, radius, angle, p, n, :TE), f, μᵣ, εᵣ)
tm_sector_fields(r, ϕ, radius, angle, p, n, f, μᵣ, εᵣ) = sector_fields(:TM, r, ϕ, kc_sector(radius, angle, p, n, :TM), sector_modal_f(r, ϕ, radius, angle, p, n, :TM), f, μᵣ, εᵣ)
te_sector_fields(r::AbstractArray{T, N}, ϕ::AbstractArray{T, N}, radius, angle, p, n, f, μᵣ, εᵣ) where {T, N} = te_sector_fields.(r, ϕ, radius, angle, p, n, f, μᵣ, εᵣ)
tm_sector_fields(r::AbstractArray{T, N}, ϕ::AbstractArray{T, N}, radius, angle, p, n, f, μᵣ, εᵣ) where {T, N} = tm_sector_fields.(r, ϕ, radius, angle, p, n, f, μᵣ, εᵣ)
te_annular_sector_fields(r, ϕ, inner_radius, outer_radius, angle, p, n, f, μᵣ, εᵣ) = sector_fields(:TE, r, ϕ, kc_annular_sector(inner_radius, outer_radius, angle, p, n, :TE), annular_sector_modal_f(r, ϕ, inner_radius, outer_radius, angle, p, n, :TE), f, μᵣ, εᵣ)
tm_annular_sector_fields(r, ϕ, inner_radius, outer_radius, angle, p, n, f, μᵣ, εᵣ) = sector_fields(:TM, r, ϕ, kc_annular_sector(inner_radius, outer_radius, angle, p, n, :TM), annular_sector_modal_f(r, ϕ, inner_radius, outer_radius, angle, p, n, :TM), f, μᵣ, εᵣ)
te_annular_sector_fields(r::AbstractArray{T, N}, ϕ::AbstractArray{T, N}, inner_radius, outer_radius, angle, p, n, f, μᵣ, εᵣ) where {T, N} = te_annular_sector_fields.(r, ϕ, inner_radius, outer_radius, angle, p, n, f, μᵣ, εᵣ)
tm_annular_sector_fields(r::AbstractArray{T, N}, ϕ::AbstractArray{T, N}, inner_radius, outer_radius, angle, p, n, f, μᵣ, εᵣ) where {T, N} = tm_annular_sector_fields.(r, ϕ, inner_radius, outer_radius, angle, p, n, f, μᵣ, εᵣ)

function angular_integral_sector(angle, p, kind)
    kind == :TE && p == 0 && return angle
    return angle / 2
end

function bessel_radial_lommel_integral(ν, x, Am, Bm)
    return Am^2 * lommel_integral(besselj, besselj, ν, x) + Bm^2 * lommel_integral(bessely, bessely, ν, x) + 2 * Am * Bm * lommel_integral(besselj, bessely, ν, x)
end

function bessel_radial_boundary_term(ν, x, Am, Bm)
    R = Am * besselj(ν, x) + Bm * bessely(ν, x)
    Rx = Am * besselj_prime(ν, x) + Bm * bessely_prime(ν, x)
    return x * R * Rx
end

function radial_energy_sector(radius, angle, p, n, kind)
    kc = kc_sector(radius, angle, p, n, kind)
    ν = sector_order(p, angle)
    x = kc * radius
    return bessel_radial_boundary_term(ν, x, 1.0, 0.0) + lommel_integral(besselj, besselj, ν, x)
end

function radial_energy_annular_sector(inner_radius, outer_radius, angle, p, n, kind)
    kc = kc_annular_sector(inner_radius, outer_radius, angle, p, n, kind)
    ν = sector_order(p, angle)
    Am, Bm = kind == :TE ? annular_sector_boundary_coeff_te(ν, inner_radius, kc) : annular_sector_boundary_coeff_tm(ν, inner_radius, kc)
    xa = kc * outer_radius
    xb = kc * inner_radius
    boundary = bessel_radial_boundary_term(ν, xa, Am, Bm) - bessel_radial_boundary_term(ν, xb, Am, Bm)
    return boundary + bessel_radial_lommel_integral(ν, xa, Am, Bm) - bessel_radial_lommel_integral(ν, xb, Am, Bm)
end

function sector_normalization(kind, kc, angle, p, radial_energy, β, f, μᵣ, εᵣ)
    ω = 2 * π * f
    factor = kind == :TE ? ω * μᵣ * _μₒ * β / kc^4 : ω * εᵣ * _εₒ * β / kc^4
    power = 0.5 * factor * angular_integral_sector(angle, p, kind) * radial_energy
    return sqrt(1 / power)
end

function te_normalization_sector(radius, angle, p, n, kc, β, f, μᵣ, εᵣ)
    return sector_normalization(:TE, kc, angle, p, radial_energy_sector(radius, angle, p, n, :TE), β, f, μᵣ, εᵣ)
end

function tm_normalization_sector(radius, angle, p, n, kc, β, f, μᵣ, εᵣ)
    return sector_normalization(:TM, kc, angle, p, radial_energy_sector(radius, angle, p, n, :TM), β, f, μᵣ, εᵣ)
end

function te_normalization_annular_sector(inner_radius, outer_radius, angle, p, n, kc, β, f, μᵣ, εᵣ)
    return sector_normalization(:TE, kc, angle, p, radial_energy_annular_sector(inner_radius, outer_radius, angle, p, n, :TE), β, f, μᵣ, εᵣ)
end

function tm_normalization_annular_sector(inner_radius, outer_radius, angle, p, n, kc, β, f, μᵣ, εᵣ)
    return sector_normalization(:TM, kc, angle, p, radial_energy_annular_sector(inner_radius, outer_radius, angle, p, n, :TM), β, f, μᵣ, εᵣ)
end

function first_n_modes_sector(N, radius, angle)
    N <= 0 && return Tuple{Symbol, Int, Int, Float64}[]
    modes = Tuple{Symbol, Int, Int, Float64}[]
    pmax = max(2, ceil(Int, sqrt(N)) + 2)
    nmax = max(2, ceil(Int, sqrt(N)) + 2)
    while true
        empty!(modes)
        for p in 0:pmax, n in 1:nmax
            push!(modes, (:TE, p, n, kc_sector(radius, angle, p, n, :TE)))
            p > 0 && push!(modes, (:TM, p, n, kc_sector(radius, angle, p, n, :TM)))
        end
        sort!(modes, by=x -> (x[4], x[1] == :TE ? 0 : 1, x[2], x[3]))
        if length(modes) >= N
            selected = modes[1:N]
            touches_pmax = any(mode -> mode[2] == pmax, selected)
            touches_nmax = any(mode -> mode[3] == nmax, selected)
            !(touches_pmax || touches_nmax) && break
            touches_pmax && (pmax *= 2)
            touches_nmax && (nmax *= 2)
        else
            pmax *= 2
            nmax *= 2
        end
    end
    return modes[1:N]
end

function first_n_modes_annular_sector(N, inner_radius, outer_radius, angle)
    N <= 0 && return Tuple{Symbol, Int, Int, Float64}[]
    modes = Tuple{Symbol, Int, Int, Float64}[]
    pmax = max(2, ceil(Int, sqrt(N)) + 2)
    nmax = max(2, ceil(Int, sqrt(N)) + 2)
    while true
        empty!(modes)
        for p in 0:pmax, n in 1:nmax
            push!(modes, (:TE, p, n, kc_annular_sector(inner_radius, outer_radius, angle, p, n, :TE)))
            p > 0 && push!(modes, (:TM, p, n, kc_annular_sector(inner_radius, outer_radius, angle, p, n, :TM)))
        end
        sort!(modes, by=x -> (x[4], x[1] == :TE ? 0 : 1, x[2], x[3]))
        if length(modes) >= N
            selected = modes[1:N]
            touches_pmax = any(mode -> mode[2] == pmax, selected)
            touches_nmax = any(mode -> mode[3] == nmax, selected)
            !(touches_pmax || touches_nmax) && break
            touches_pmax && (pmax *= 2)
            touches_nmax && (nmax *= 2)
        else
            pmax *= 2
            nmax *= 2
        end
    end
    return modes[1:N]
end
