const _εₒ::Float64 = 8.8541878188e−12
const _μₒ::Float64 = 1.25663706127e−6 
const _c::Float64 = 1/sqrt(_μₒ * _εₒ)
const _η::Float64 = sqrt(_μₒ/_εₒ)

"""
    wavenumber(f, μᵣ, εᵣ)

Computes the wavenumber ω√(μ⋅ε).

# Arguments
- `f`: Frequency in Hz
- `μᵣ`: Relative permeability
- `εᵣ`: Relative permittivity
"""
function wavenumber(f, μᵣ, εᵣ)
    μ = μᵣ * _μₒ
    ε = εᵣ * _εₒ
    return 2π * f * sqrt(μ * ε)
end

"""
    cutoff_frequency(kc, μᵣ, εᵣ)

Computes the cutoff frequency of a mode given the cutoff wavenumber `kc`.

# Arguments
- `kc`: Cutoff wavenumber
- `μᵣ`: Relative permeability
- `εᵣ`: Relative permittivity
"""
function cutoff_frequency(kc, μᵣ, εᵣ)
    μ = μᵣ * _μₒ
    ε = εᵣ * _εₒ
    return 1/(2π * sqrt(μ * ε)) * kc
end


"""
    propagation_constant(kc, f, μᵣ, εᵣ)

Computes γ = α + jβ given the cutoff wavenumber `kc`.

# Arguments
- `kc`: Cutoff wavenumber
- `f`: Frequency in Hz
- `μᵣ`: Relative permeability
- `εᵣ`: Relative permittivity
"""
function propagation_constant(kc, f, μᵣ, εᵣ)
    μ = μᵣ * _μₒ
    ε = εᵣ * _εₒ
    ω = 2π * f
    k = ω * sqrt(μ * ε)
    return im * sqrt(Complex(kc^2 - k^2)) # general form
end

"""
    phase_constant(kc, f, μᵣ, εᵣ)

Computes β given the cutoff wavenumber `kc`.

# Arguments
- `kc`: Cutoff wavenumber
- `f`: Frequency in Hz
- `μᵣ`: Relative permeability
- `εᵣ`: Relative permittivity
"""
function phase_constant(kc, f, μᵣ, εᵣ)
    μ = μᵣ * _μₒ
    ε = εᵣ * _εₒ
    ω = 2π * f
    k = ω * sqrt(μ * ε)
    return sqrt(Complex(k^2 - kc^2))
end

"""
    attenuation_factor(kc, f, μᵣ, εᵣ)

Computes α, the attenuation constant of a given mode for frequencies below the cutoff.

# Arguments
- `kc`: Cutoff wavenumber
- `f`: Frequency in Hz
- `μᵣ`: Relative permeability
- `εᵣ`: Relative permittivity
"""
function attenuation_factor(kc, f, μᵣ, εᵣ)
    return real(propagation_constant(kc, f, μᵣ, εᵣ))
end


"""
    mode_wavelength(β)

Returns the guided wavelength (mode wavelength) of the propagating mode given the phase constant `β`.

# Arguments
- `β`: Phase constant
"""
function mode_wavelength(β)
    return ifelse(imag(β) == 0.0, 2π/real(β), 0.0)
end


"""
    te_coefficients(kc, β, f, μᵣ, εᵣ) 

Calculates the transverse electric field (E) and magnetic field (H) coupling coefficients
for a **Transverse Electric (TE)** mode.
"""
function te_coefficients(kc, β, f, μᵣ, εᵣ)
    μ = μᵣ * _μₒ
    ω = 2 * π * f
    return (Complex(ω*μ/kc^2), Complex(β/kc^2))
end


"""
    tm_coefficients(kc, β, f, μᵣ, εᵣ)

Calculates the transverse electric field (E) and magnetic field (H) coupling coefficients
for a **Transverse Magnetic (TM)** mode.
"""
function tm_coefficients(kc, β, f, μᵣ, εᵣ)
    ε = εᵣ * _εₒ 
    ω = 2 * π * f
    return (Complex(β / kc^2), Complex(ω * ε / kc^2))
end

# Impedance

"""
    mode_impedance(T, fcutoff, f, μᵣ, εᵣ)

Computes the mode impedance of a T ∈ (:TE, :TM) mode given the cutoff frequency.

# Arguments
- `T`: Mode type (`:TE` or `:TM`)
- `fcutoff`: Cutoff frequency
- `f`: Frequency in Hz
- `μᵣ`: Relative permeability
- `εᵣ`: Relative permittivity
"""
function mode_impedance(T, fcutoff, f, μᵣ, εᵣ)
    if T == :TE
        return _η / sqrt(μᵣ * εᵣ) / sqrt(Complex(1 - (fcutoff / f)^2))
    elseif T == :TM
        return _η / sqrt(μᵣ * εᵣ) * sqrt(Complex(1 - (fcutoff / f)^2))
    else
        @error "$T not suported."
    end
end

# Derivatives of Bessel functions

"""
    besselj_prime(m, x)

Bessel J function first derivative.
"""
besselj_prime(m, x) = (besselj(m - 1, x) - besselj(m + 1, x)) / 2

"""
    bessely_prime(m, x)

Bessel Y function first derivative.
"""
bessely_prime(m, x) = (bessely(m - 1, x) - bessely(m + 1, x)) / 2

"""
    hankelh1_prime(m, x)

Hankel H¹ function first derivative.
"""
hankelh1_prime(m, x) = (hankelh1(m - 1, x) - hankelh1(m + 1, x)) / 2

"""
    hankelh1_prime(m, x)

Hankel H² function first derivative.
"""
hankelh2_prime(m, x) = (hankelh2(m - 1, x) - hankelh2(m + 1, x)) / 2


"""
    besselj_doubleprime(m, x)

Second derivative of bessel J function of order m.
"""
function besselj_doubleprime(m, x)
    return -(1 / x) * besselj_prime(m, x) - (1 - (m^2) / (x^2)) * besselj(m, x)
end

"""
    bessely_doubleprime(m, x)

Second derivative of bessel Y function of order m.
"""
function bessely_doubleprime(m, x)
    return -(1 / x) * bessely_prime(m, x) - (1 - (m^2) / (x^2)) * bessely(m, x)
end


# Integrals related to Bessel functions

"""
    lommel_integral(C::F1, D::F2, m, x) where {F1, F2} -> Real

Calculates the integral of the product of two Bessel functions of the same order
and argument including
"""
function lommel_integral(C::F1, D::F2, m, x) where {F1, F2}
    return x^2/4 * (2*C(m, x) * D(m, x) - C(m-1, x) * D(m+1, x) - C(m+1, x) * D(m-1, x))
end

"""
    lommel_integral_frac(C::F1, D::F2, m, x) where {F1, F2} -> Real

Calculates the integral of the product of two Bessel functions of the same order
and argument including the factor 1/x.
"""
function lommel_integral_frac(C::F1, D::F2, m, x) where {F1, F2}
    return 1/(2*m) * (C(m, x)*D(m-1, x) - C(m-1, x)*D(m, x))
end

"""
    lommel_integral_frac_u_v(C::F1, D::F2, m, x) where {F1, F2} -> Real

Calculates the integral of the product of two Bessel functions of order u and v including the factor 1/x.
"""
function lommel_integral_frac_u_v(C::F1, D::F2, u, v, x) where {F1, F2}
    return -(x * (C(u+1, x)*D(v, x) - C(u, x) * D(v+1, x) ))/(u^2-v^2) + (C(u,x)*D(v,x))/(u+v)
end

"""
    lommel_integral_uv(C::F1, D::F2, m, x) where {F1, F2}

Calculates the integral of the product of two Bessel functions of order u and v.
"""
function lommel_integral_uv(C::F1, D::F2, u, v, x) where {F1, F2}
    return x / (u^2-v^2) * (u*C(u+1,x)*D(v, x) - v*C(u, x)*D(v+1, x))
end

