"""
    spheroidal_scale_factors(a, ξ, η; oblate=false)

Return spheroidal metric factors `(hξ, hη, hϕ)` for prolate (`oblate=false`) or
oblate (`oblate=true`) coordinates with semifocal distance `a`.
"""
function spheroidal_scale_factors(a, ξ, η; oblate::Bool = false)
    if oblate
        return scale_factors_oblate(a, ξ, η)
    end
    return scale_factors_prolate(a, ξ, η)
end

@inline function _trapz_periodic(vals, Δ)
    return Δ * sum(vals)
end

"""
    spheroidal_power_density_xi(Eη, Eϕ, Hη, Hϕ)

Time-averaged Poynting flux component along `êξ`:

`Sξ = 0.5*real(Eη*conj(Hϕ) - Eϕ*conj(Hη))`
"""
@inline function spheroidal_power_density_xi(Eη, Eϕ, Hη, Hϕ)
    return 0.5 * real(Eη * conj(Hϕ) - Eϕ * conj(Hη))
end

"""
    spheroidal_power_on_xi_surface(Efun, Hfun, a, ξ0; oblate=false, Nη=256, Nϕ=256)

Numerically integrate the net time-averaged power crossing the closed surface `ξ=ξ0`.

`Efun(η, ϕ)` and `Hfun(η, ϕ)` must return spheroidal components:
`(Eξ, Eη, Eϕ)` and `(Hξ, Hη, Hϕ)`.

Integration is performed with periodic trapezoidal rules in `η ∈ [-1, 1]` and
`ϕ ∈ [0, 2π)`. Surface element on `ξ=const` is `dS = hη*hϕ dη dϕ`.
"""
function spheroidal_power_on_xi_surface(Efun, Hfun, a, ξ0;
                                        oblate::Bool = false,
                                        Nη::Int = 256,
                                        Nϕ::Int = 256)
    Nη < 2 && throw(ArgumentError("Nη must be >= 2"))
    Nϕ < 2 && throw(ArgumentError("Nϕ must be >= 2"))

    ηs = range(-1.0, 1.0; length = Nη)
    ϕs = range(0.0, 2π; length = Nϕ + 1)[1:end-1]
    Δη = step(ηs)
    Δϕ = 2π / Nϕ

    total = 0.0
    for η in ηs
        _, hη, hϕ = spheroidal_scale_factors(a, ξ0, η; oblate = oblate)
        row = 0.0
        for ϕ in ϕs
            _, Eη, Eϕ = Efun(η, ϕ)
            _, Hη, Hϕ = Hfun(η, ϕ)
            Sξ = spheroidal_power_density_xi(Eη, Eϕ, Hη, Hϕ)
            row += Sξ * hη * hϕ
        end
        total += row * Δϕ
    end
    return total * Δη
end

@inline _is_oblate_mode(::SpheroidalB{I, T, T}) where {I, T} = false
@inline _is_oblate_mode(::SpheroidalB{I, T, Complex{T}}) where {I, T} = true

"""
    spheroidal_mn_power_on_xi_surface(mode::SpheroidalB, k, ξ0;
                                      family=:z, even=true, oblate=_is_oblate_mode(mode),
                                      radial=4, CE=1+0im, CH=1+0im,
                                      E_from=:M, H_from=:N, Nη=256, Nϕ=256)

Power flux on `ξ=ξ0` using hardcoded spheroidal vectors `M/N`.

Field model:
- If `E_from=:M`, `E = CE*M`; if `E_from=:N`, `E = CE*N`.
- If `H_from=:M`, `H = CH*M`; if `H_from=:N`, `H = CH*N`.

This keeps the TE/TM convention explicit via `CE`, `CH`, and source selectors.
"""
function spheroidal_mn_power_on_xi_surface(mode::SpheroidalB, k, ξ0;
                                           family::Symbol = :z,
                                           even::Bool = true,
                                           oblate::Bool = _is_oblate_mode(mode),
                                           radial::Int = 4,
                                           CE = 1.0 + 0.0im,
                                           CH = 1.0 + 0.0im,
                                           E_from::Symbol = :M,
                                           H_from::Symbol = :N,
                                           Nη::Int = 256,
                                           Nϕ::Int = 256)
    a = abs(mode.c / k)

    Efun = function (η, ϕ)
        Mξ, Mη, Mϕ, Nξ, Nη, Nϕ = mn_spheroidal_vector(ξ0, η, ϕ, mode, k;
                                                      family = family, even = even,
                                                      oblate = oblate, radial = radial)
        if E_from == :M
            return (CE * Mξ, CE * Mη, CE * Mϕ)
        elseif E_from == :N
            return (CE * Nξ, CE * Nη, CE * Nϕ)
        end
        throw(ArgumentError("E_from must be :M or :N"))
    end

    Hfun = function (η, ϕ)
        Mξ, Mη, Mϕ, Nξ, Nη, Nϕ = mn_spheroidal_vector(ξ0, η, ϕ, mode, k;
                                                      family = family, even = even,
                                                      oblate = oblate, radial = radial)
        if H_from == :M
            return (CH * Mξ, CH * Mη, CH * Mϕ)
        elseif H_from == :N
            return (CH * Nξ, CH * Nη, CH * Nϕ)
        end
        throw(ArgumentError("H_from must be :M or :N"))
    end

    return spheroidal_power_on_xi_surface(Efun, Hfun, a, ξ0;
                                          oblate = oblate, Nη = Nη, Nϕ = Nϕ)
end

"""
    spheroidal_mn_normalization_on_xi_surface(mode::SpheroidalB, k, ξ0;
                                              target=1.0, kwargs...)

Return multiplicative normalization `A` such that `|A|^2 * P = target`, where `P`
is computed by [`spheroidal_mn_power_on_xi_surface`](@ref).
"""
function spheroidal_mn_normalization_on_xi_surface(mode::SpheroidalB, k, ξ0;
                                                   target::Real = 1.0, kwargs...)
    P = spheroidal_mn_power_on_xi_surface(mode, k, ξ0; kwargs...)
    abs(P) == 0 && throw(ArgumentError("Computed power is zero; cannot normalize."))
    return sqrt(target / abs(P))
end

@inline function _simpson_uniform(vals, h)
    n = length(vals)
    n < 3 && throw(ArgumentError("Need at least 3 samples for Simpson integration"))
    iseven(n) && throw(ArgumentError("Simpson requires odd number of samples"))
    s = vals[1] + vals[end]
    s += 4 * sum(vals[2:2:end-1])
    s += 2 * sum(vals[3:2:end-2])
    return (h / 3) * s
end

@inline function _recommended_nphi(m::Int; safety::Int = 4)
    # Max harmonic from x/y families is roughly m±1; z/r are lower.
    n = 2 * (m + 1 + safety)
    return max(32, n)
end

"""
    spheroidal_power_eta_kernel(mode::SpheroidalB, k, ξ0, η;
                                family=:z, even=true, oblate=_is_oblate_mode(mode),
                                radial=4, CE=1+0im, CH=1+0im, E_from=:M, H_from=:N,
                                Nϕ=nothing)

Compute the 1D kernel

`K(η) = ∫₀²π Sξ(ξ0,η,ϕ) hη(ξ0,η) hϕ(ξ0,η) dϕ`

using periodic trapezoidal quadrature in `ϕ`. For trigonometric-polynomial dependence
of the hardcoded spheroidal families, this is spectrally accurate and effectively exact
once `Nϕ` resolves the highest harmonic.
"""
function spheroidal_power_eta_kernel(mode::SpheroidalB, k, ξ0, η;
                                     family::Symbol = :z,
                                     even::Bool = true,
                                     oblate::Bool = _is_oblate_mode(mode),
                                     radial::Int = 4,
                                     CE = 1.0 + 0.0im,
                                     CH = 1.0 + 0.0im,
                                     E_from::Symbol = :M,
                                     H_from::Symbol = :N,
                                     Nϕ = nothing)
    a = abs(mode.c / k)
    _, hη, hϕ = spheroidal_scale_factors(a, ξ0, η; oblate = oblate)
    nphi = isnothing(Nϕ) ? _recommended_nphi(mode.m) : Int(Nϕ)
    nphi < 4 && throw(ArgumentError("Nϕ must be >= 4"))

    ϕs = range(0.0, 2π; length = nphi + 1)[1:end-1]
    Δϕ = 2π / nphi
    acc = 0.0
    for ϕ in ϕs
        Mξ, Mη, Mϕ, Nξ, Nη, Nϕ = mn_spheroidal_vector(ξ0, η, ϕ, mode, k;
                                                      family = family, even = even,
                                                      oblate = oblate, radial = radial)
        Eξ, Eη, Eϕ = if E_from == :M
            (CE * Mξ, CE * Mη, CE * Mϕ)
        elseif E_from == :N
            (CE * Nξ, CE * Nη, CE * Nϕ)
        else
            throw(ArgumentError("E_from must be :M or :N"))
        end
        Hξ, Hη, Hϕ = if H_from == :M
            (CH * Mξ, CH * Mη, CH * Mϕ)
        elseif H_from == :N
            (CH * Nξ, CH * Nη, CH * Nϕ)
        else
            throw(ArgumentError("H_from must be :M or :N"))
        end
        acc += spheroidal_power_density_xi(Eη, Eϕ, Hη, Hϕ)
    end
    return acc * Δϕ * hη * hϕ
end

"""
    spheroidal_mn_power_on_xi_surface_1d(mode::SpheroidalB, k, ξ0;
                                         family=:z, even=true, oblate=_is_oblate_mode(mode),
                                         radial=4, CE=1+0im, CH=1+0im,
                                         E_from=:M, H_from=:N,
                                         Nη=257, Nϕ=nothing)

Semi-analytical (quasi-analytical) modal power on `ξ=ξ0`:
- `ϕ` integral is performed via periodic spectral trapezoidal rule (exact for resolved
  finite trigonometric content),
- remaining integral is 1D in `η`, integrated with composite Simpson.
"""
function spheroidal_mn_power_on_xi_surface_1d(mode::SpheroidalB, k, ξ0;
                                              family::Symbol = :z,
                                              even::Bool = true,
                                              oblate::Bool = _is_oblate_mode(mode),
                                              radial::Int = 4,
                                              CE = 1.0 + 0.0im,
                                              CH = 1.0 + 0.0im,
                                              E_from::Symbol = :M,
                                              H_from::Symbol = :N,
                                              Nη::Int = 257,
                                              Nϕ = nothing)
    Nη < 3 && throw(ArgumentError("Nη must be >= 3"))
    iseven(Nη) && throw(ArgumentError("Nη must be odd for Simpson integration"))
    ηs = range(-1.0, 1.0; length = Nη)
    Δη = step(ηs)

    K = similar(collect(ηs))
    for i in eachindex(ηs)
        η = ηs[i]
        K[i] = spheroidal_power_eta_kernel(mode, k, ξ0, η;
                                           family = family, even = even, oblate = oblate,
                                           radial = radial, CE = CE, CH = CH,
                                           E_from = E_from, H_from = H_from, Nϕ = Nϕ)
    end
    return _simpson_uniform(K, Δη)
end

"""
    spheroidal_mn_normalization_on_xi_surface_1d(mode::SpheroidalB, k, ξ0;
                                                  target=1.0, kwargs...)

Normalization factor based on [`spheroidal_mn_power_on_xi_surface_1d`](@ref):
returns `A` such that `|A|^2*P = target`.
"""
function spheroidal_mn_normalization_on_xi_surface_1d(mode::SpheroidalB, k, ξ0;
                                                      target::Real = 1.0, kwargs...)
    P = spheroidal_mn_power_on_xi_surface_1d(mode, k, ξ0; kwargs...)
    abs(P) == 0 && throw(ArgumentError("Computed power is zero; cannot normalize."))
    return sqrt(target / abs(P))
end

"""
    spheroidal_power_eta_kernels(mode::SpheroidalB, k, ξ0, η;
                                 family=:z, even=true, oblate=_is_oblate_mode(mode),
                                 radial=4, Nϕ=nothing)

Compute complex ϕ-integrated kernels at fixed `(ξ0, η)`:

- `JMM = ∫ (Mη*conj(Mϕ) - Mϕ*conj(Mη)) hηhϕ dϕ`
- `JMN = ∫ (Mη*conj(Nϕ) - Mϕ*conj(Nη)) hηhϕ dϕ`
- `JNM = ∫ (Nη*conj(Mϕ) - Nϕ*conj(Mη)) hηhϕ dϕ`
- `JNN = ∫ (Nη*conj(Nϕ) - Nϕ*conj(Nη)) hηhϕ dϕ`

These are the spheroidal analog of reusable angular kernels: once computed, power for
any linear combination of `M/N` with complex coefficients is reconstructed without
re-integrating in `ϕ`.
"""
function spheroidal_power_eta_kernels(mode::SpheroidalB, k, ξ0, η;
                                      family::Symbol = :z,
                                      even::Bool = true,
                                      oblate::Bool = _is_oblate_mode(mode),
                                      radial::Int = 4,
                                      Nϕ = nothing)
    a = abs(mode.c / k)
    _, hη, hϕ = spheroidal_scale_factors(a, ξ0, η; oblate = oblate)
    nphi = isnothing(Nϕ) ? _recommended_nphi(mode.m) : Int(Nϕ)
    nphi < 4 && throw(ArgumentError("Nϕ must be >= 4"))
    ϕs = range(0.0, 2π; length = nphi + 1)[1:end-1]
    Δϕ = 2π / nphi

    JMM = 0.0 + 0.0im
    JMN = 0.0 + 0.0im
    JNM = 0.0 + 0.0im
    JNN = 0.0 + 0.0im

    for ϕ in ϕs
        Mξ, Mη, Mϕ, Nξ, Nη, Nϕ = mn_spheroidal_vector(ξ0, η, ϕ, mode, k;
                                                      family = family, even = even,
                                                      oblate = oblate, radial = radial)
        JMM += (Mη * conj(Mϕ) - Mϕ * conj(Mη))
        JMN += (Mη * conj(Nϕ) - Mϕ * conj(Nη))
        JNM += (Nη * conj(Mϕ) - Nϕ * conj(Mη))
        JNN += (Nη * conj(Nϕ) - Nϕ * conj(Nη))
    end

    scale = Δϕ * hη * hϕ
    return (JMM = JMM * scale, JMN = JMN * scale, JNM = JNM * scale, JNN = JNN * scale)
end

"""
    spheroidal_power_eta_kernels_grid(mode::SpheroidalB, k, ξ0, ηs;
                                      family=:z, even=true, oblate=_is_oblate_mode(mode),
                                      radial=4, Nϕ=nothing)

Compute reusable ϕ-integrated kernels on an η grid.
Returns `(η=ηs, JMM, JMN, JNM, JNN)`.
"""
function spheroidal_power_eta_kernels_grid(mode::SpheroidalB, k, ξ0, ηs;
                                           family::Symbol = :z,
                                           even::Bool = true,
                                           oblate::Bool = _is_oblate_mode(mode),
                                           radial::Int = 4,
                                           Nϕ = nothing)
    n = length(ηs)
    JMM = Vector{ComplexF64}(undef, n)
    JMN = Vector{ComplexF64}(undef, n)
    JNM = Vector{ComplexF64}(undef, n)
    JNN = Vector{ComplexF64}(undef, n)
    for i in eachindex(ηs)
        K = spheroidal_power_eta_kernels(mode, k, ξ0, ηs[i];
                                         family = family, even = even, oblate = oblate,
                                         radial = radial, Nϕ = Nϕ)
        JMM[i] = K.JMM
        JMN[i] = K.JMN
        JNM[i] = K.JNM
        JNN[i] = K.JNN
    end
    return (η = collect(ηs), JMM = JMM, JMN = JMN, JNM = JNM, JNN = JNN)
end

@inline function _pick_kernel(kernels, E_from::Symbol, H_from::Symbol, i::Int)
    if E_from == :M && H_from == :M
        return kernels.JMM[i]
    elseif E_from == :M && H_from == :N
        return kernels.JMN[i]
    elseif E_from == :N && H_from == :M
        return kernels.JNM[i]
    elseif E_from == :N && H_from == :N
        return kernels.JNN[i]
    end
    throw(ArgumentError("E_from/H_from must be :M or :N"))
end

"""
    spheroidal_modal_power_from_eta_kernels(kernels;
                                            CE=1+0im, CH=1+0im,
                                            E_from=:M, H_from=:N)

Reconstruct modal power from precomputed η-kernels:

`P = 0.5 ∫ Re(CE*conj(CH) * J_{E,H}(η)) dη`

with `J_{E,H} ∈ {JMM, JMN, JNM, JNN}` depending on `E_from/H_from`.
"""
function spheroidal_modal_power_from_eta_kernels(kernels;
                                                 CE = 1.0 + 0.0im,
                                                 CH = 1.0 + 0.0im,
                                                 E_from::Symbol = :M,
                                                 H_from::Symbol = :N)
    η = kernels.η
    n = length(η)
    n < 3 && throw(ArgumentError("Need at least 3 η nodes"))
    iseven(n) && throw(ArgumentError("η grid length must be odd (Simpson)"))
    Δη = η[2] - η[1]
    pref = CE * conj(CH)

    vals = Vector{Float64}(undef, n)
    for i in eachindex(η)
        J = _pick_kernel(kernels, E_from, H_from, i)
        vals[i] = 0.5 * real(pref * J)
    end
    return _simpson_uniform(vals, Δη)
end

"""
    spheroidal_mn_power_on_xi_surface_1d_from_kernels(mode::SpheroidalB, k, ξ0;
                                                       family=:z, even=true,
                                                       oblate=_is_oblate_mode(mode),
                                                       radial=4, Nη=257, Nϕ=nothing,
                                                       CE=1+0im, CH=1+0im,
                                                       E_from=:M, H_from=:N)

Semi-analytical power via reusable η-kernels (spheroidal counterpart of Av/Avp-style
decomposition): first build ϕ-integrated kernels, then integrate only in η.
"""
function spheroidal_mn_power_on_xi_surface_1d_from_kernels(mode::SpheroidalB, k, ξ0;
                                                            family::Symbol = :z,
                                                            even::Bool = true,
                                                            oblate::Bool = _is_oblate_mode(mode),
                                                            radial::Int = 4,
                                                            Nη::Int = 257,
                                                            Nϕ = nothing,
                                                            CE = 1.0 + 0.0im,
                                                            CH = 1.0 + 0.0im,
                                                            E_from::Symbol = :M,
                                                            H_from::Symbol = :N)
    Nη < 3 && throw(ArgumentError("Nη must be >= 3"))
    iseven(Nη) && throw(ArgumentError("Nη must be odd for Simpson integration"))
    ηs = range(-1.0, 1.0; length = Nη)
    kernels = spheroidal_power_eta_kernels_grid(mode, k, ξ0, ηs;
                                                family = family, even = even,
                                                oblate = oblate, radial = radial, Nϕ = Nϕ)
    return spheroidal_modal_power_from_eta_kernels(kernels;
                                                   CE = CE, CH = CH,
                                                   E_from = E_from, H_from = H_from)
end
