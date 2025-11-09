include("./coax_functions.jl")

"""
    characteristic_coax_equation_te(m, kc, a, b)

The zeros of this function correspond to the cutoff wave number of the TE modes.

# Arguments
- `m`: Azimuthal mode index
- `kc`: Cutoff wavenumber
- `a`: Outer radius
- `b`: Inner radius
"""
characteristic_coax_equation_te(m, kc, a, b) = besselj_prime(m, kc * a) * bessely_prime(m, kc * b) - bessely_prime(m, kc * a) * besselj_prime(m, kc * b)
"""
    characteristic_coax_equation_tm(m, kc, a, b)

The zeros of this function correspond to the cutoff wave number of the TM modes.

# Arguments
- `m`: Azimuthal mode index
- `kc`: Cutoff wavenumber
- `a`: Outer radius
- `b`: Inner radius
"""
characteristic_coax_equation_tm(m, kc, a, b) = besselj(m, kc * a) * bessely(m, kc * b) - bessely(m, kc * a) * besselj(m, kc * b)

"""
    characteristic_coax_equation_tm_prime(m, kc, a, b)

Derivative of characteristic_coax_equation_tm(m, kc, a, b)
"""
function characteristic_coax_equation_tm_prime(m, kc, a, b)
    return a * besselj_prime(m, kc * a) * bessely(m, kc * b) + b * besselj(m, kc * a) * bessely_prime(m, kc * b) -
           b * bessely(m, kc * a) * besselj_prime(m, kc * b) - a * bessely_prime(m, kc * a) * besselj(m, kc * b)
end


"""
    characteristic_coax_equation_te_prime(m, kc, a, b)

Derivative of characteristic_coax_equation_te(m, kc, a, b)
"""
function characteristic_coax_equation_te_prime(m, kc, a, b)
    return a * besselj_doubleprime(m, kc * a) * bessely_prime(m, kc * b) + b * besselj_prime(m, kc * a) * bessely_doubleprime(m, kc * b) -
           b * bessely_prime(m, kc * a) * besselj_doubleprime(m, kc * b) - a * bessely_doubleprime(m, kc * a) * besselj_prime(m, kc * b)
end


"""
    zps(m, n, ratio)

- Cochran approximation for TM modes in coaxial waveguides.
- Cochran, J.A.. (1963). Further Formulas for Calculating Approximate Values of the Zeros of Certain Combinations of Bessel Functions. 10.1109/TMTT.1963.1125722
"""
zps(m, n, ratio) = sqrt((n * π)^2 / (ratio - 1)^2 + (4 * m^2 - 1) / (ratio + 1)^2)

"""
    z¹ps(m, n, ratio)

- Cochran approximation for TE modes in coaxial waveguides.
- Cochran, J.A.. (1963). Further Formulas for Calculating Approximate Values of the Zeros of Certain Combinations of Bessel Functions. 10.1109/TMTT.1963.1125722
"""
z¹ps(m, n, ratio) = ifelse(n == 0, z¹ps_0(m, n, ratio), z¹ps_n(m, n, ratio))
z¹ps_n(m, n, ratio) = sqrt((n * π)^2 / (ratio - 1)^2 + (4 * m^2 + 3) / (ratio + 1)^2)
z¹ps_0(m, n, ratio) = (2*m)/(ratio+1) * ((1+(ratio-1)^2)/(6*(ratio+1)^2))


"""
    asymtotic_coax_zero_tm(m, n, a, b)

Compute the asymptotic approximation of the nth zero of a mode TM_{m, n}.

Cross-Product

```
J(m,x)*Y(m, λ*x) - Y(m,x)*J(m, λ*x)
```

nth zero is equal to

```
α + p/α + (q-p^2)/(α^3) + (r - 4*p*q + 2*p^3)/(α^5)
```

```
α = n * π / (λ - 1)
p = (μ-1)/(8*λ)
q = ((μ-1)*(μ-25)*(λ^3-1))/(6*(4*λ)^3*(λ-1))
r = ((μ-1)*(μ^2-114*μ+1073)*(λ^5-1))/(5*(4*λ)^5*(λ-1))
```
Reference: Reference: https://dlmf.nist.gov/10.21#ii

"""
function asymtotic_coax_zero_tm(m, n, a, b)
    λ = b/a
    μ = 4 * m^2
    α = n * π / (λ - 1)
    p = (μ - 1) / (8 * λ)
    q = ((μ - 1) * (μ - 25) * (λ^3 - 1)) / (6 * (4 * λ)^3 * (λ - 1))
    r = ((μ - 1) * (μ^2 - 114 * μ + 1073) * (λ^5 - 1)) / (5 * (4 * λ)^5 * (λ - 1))
    return (α + p / α + (q - p^2) / (α^3) + (r - 4 * p * q + 2 * p^3) / (α^5))
end

"""
    asymtotic_coax_zero_te(m, n, a, b)

Compute the asymptotic approximation of the nth zero of a mode TE_{m, n}.


Cross-Product:

```
J'(m,x)*Y'(m, λ*x) - Y'(m,x)*J'(m, λ*x)
```

nth zero is equal to

```
α + p/α + (q-p^2)/(α^3) + (r - 4*p*q + 2*p^3)/(α^5)
```

where
```
α = (n-1)*π/(λ-1)
p = (μ+3)/(8*λ)
q = ((μ^2+46*μ-63)*(λ^3-1))/(6*(4*λ)^3*(λ-1))
r = (μ^3 + 185*μ^2 - 2053*μ + 1899)*(λ^5-1))/(5*(4*λ)^5*(λ-1))
```
Reference: https://dlmf.nist.gov/10.21#ii

"""
function asymtotic_coax_zero_te(m, n, a, b)
    λ = b / a
    μ = 4 * m^2
    n == 1 && return 0.0
    α = (n - 1) * π / (λ - 1)
    p = (μ + 3) / (8 * λ)
    q = ((μ^2 + 46 * μ - 63) * (λ^3 - 1)) / (6 * (4 * λ)^3 * (λ - 1))
    r = (μ^3 + 185 * μ^2 - 2053 * μ + 1899) * (λ^5 - 1) / (5 * (4 * λ)^5 * (λ - 1))
    return (α + p / α + (q - p^2) / (α^3) + (r - 4 * p * q + 2 * p^3) / (α^5)) #/ b
end

"""
    kc_coax(a, b, m, n, T)

Compute the cutoff wave number for a coaxial waveguide.

# Arguments
- `a`: Outer radius
- `b`: Inner radius
- `m`: Azimuthal mode index
- `n`: Radial mode index
- `T`: Mode type (`:TE`, `:TM`, or `:TEM`)

# Note
For TEM mode, returns 0.0 (no cutoff).
"""
function kc_coax(a, b, m, n, T)
    maxv = max(a, b)
    minv = min(a,b)
    k = maxv/minv
	if T == :TE
		return kc_coax_te(a, b, m, n) * k / maxv
	elseif T == :TM
		return kc_coax_tm(a, b, m, n) * k / maxv
	elseif T == :TEM
		return 0.0
	else
		@error "Not supported type"
	end
end


"""
    coax_boundary_coeff_te(m, b, kc)

Returns `Am` and `Bm` used on `Φ¹ = kc * (Am * besselj_prime(m, kc * r) + Bm * bessely_prime(m, kc * r))`.

# Arguments
- `m`: Azimuthal mode index
- `b`: Inner radius
- `kc`: Cutoff wavenumber

# Returns
`(Am, Bm)` where `Am = 1.0` by convention
"""
function coax_boundary_coeff_te(m, b, kc)
    Am = 1.0
    Bm = -besselj_prime(m, kc * b) / bessely_prime(m, kc * b)
    return (Am, Bm)
end

"""
    coax_boundary_coeff_tm(m, b, kc)

Returns `Am` and `Bm` used on `Φ = Am * besselj(m, kc * r) + Bm * bessely(m, kc * r)`.

# Arguments
- `m`: Azimuthal mode index
- `b`: Inner radius
- `kc`: Cutoff wavenumber

# Returns
`(Am, Bm)` where `Am = 1.0` by convention
"""
function coax_boundary_coeff_tm(m, b, kc)
    Am = 1.0
    Bm = -besselj(m, kc * b) / bessely(m, kc * b)
    return (Am, Bm)
end


"""
    coax_modal_f(r, θ, kc, Am, Bm, m)

Modal functions for a TE/TM mode on a coaxial waveguide.

# Arguments
- `r`: Radial coordinate
- `θ`: Angular coordinate
- `kc`: Cutoff wavenumber
- `Am`, `Bm`: Boundary coefficients (from `coax_boundary_coeff_te/tm`)
- `m`: Azimuthal mode index

# Returns
A tuple (∂Hz/∂ρ, ∂Hz/∂θ, Hz) for a given (r, θ).
"""
function coax_modal_f(r, θ, kc, Am, Bm, m)
    Φ = Am * besselj(m, kc * r) + Bm * bessely(m, kc * r)
    Φ¹ = kc * (Am * besselj_prime(m, kc * r) + Bm * bessely_prime(m, kc * r))

    smθ, cmθ = sincos(m * θ)

    ∂ψᵢ = Φ¹ * cmθ # ∂Fz/∂r
    ∂ψⱼ = m / r * Φ * -smθ # ∂Fz/∂θ
    ψₖ = Φ * cmθ # Fz

    return (∂ψᵢ, ∂ψⱼ, ψₖ)
end

#TODO: Vertical / Horizontal polarization (-pi/2)
"""
    tm_coax_fields(r, θ, Am, Bm, m, kc, c_e, c_h)

# Arguments

- `r`: r polar coordinate.
- `θ`: angle in polar coordinates.
- `Am`: First value returned by `coax_boundary_coeff_te`. Is always 1.0
- `Bm`: Second value retuened by `coax_boundary_coeff_te`.
- `m`:  The azimuthal mode number.
- `n`: radial mode number.
- `kc`: cutoff wavenumber.
- `c_e`: im*ω*μ/kc^2
- `c_h`: -im*β/kc^2
- `a`: outer radius.
- `b`: inner radius
- `f`: frequency.
- `μᵣ`: relative permeability
- `εᵣ`: relative permittivity

# Returns

- `(Er, Eθ, Ez, Hr, Hθ, Hz)`

Calculate the electric and magnetic field components in polar coordinates of a TM_{m,n} mode given the dimensions of the circular section with radius r. 
It is assumed that the center of the circular section is at (0,0).

"""
function tm_coax_fields(r, θ, Am, Bm, m, kc, c_e, c_h)

    ∂ψᵢ, ∂ψⱼ, ψₖ = coax_modal_f(r, θ, kc, Am, Bm, m)

    Er = -c_e * ∂ψᵢ
    Eθ = -c_e * ∂ψⱼ
    Ez =  -im *  ψₖ
    Hr = +c_h * ∂ψⱼ
    Hθ = -c_h * ∂ψᵢ
    Hz = zero(Ez)

    return (Er, Eθ, Ez, Hr, Hθ, Hz)
end


"""
    tm_coax_fields(r, θ, a, b, m, n, f, μᵣ, εᵣ)

# Arguments

- `r`: r polar coordinate.
- `θ`: angle in polar coordinates.
- `a`: outer radius.
- `b`: inner radius
- `m`: azimuthal mode number.
- `n`: radial mode number.
- `f`: frequency.
- `μᵣ`: relative permeability
- `εᵣ`: relative permittivity

# Returns

- `(Er, Eθ, Ez, Hr, Hθ, Hz)`

Calculate the electric and magnetic field components in polar coordinates of a TM_{m,n} mode given the dimensions of the circular section with radius r. 
It is assumed that the center of the circular section is at (0,0).

"""
function tm_coax_fields(r, θ, a, b, m, n, f, μᵣ, εᵣ)

    kc = kc_coax(a, b, m, n, :TM)
    β = phase_constant(kc, f, μᵣ, εᵣ)
    Am, Bm = coax_boundary_coeff_tm(m, b, kc)
    c_e, c_h = tm_coefficients(kc, β, f, μᵣ, εᵣ)

    (Er, Eθ, Ez, Hr, Hθ, Hz) = tm_coax_fields(r, θ, Am, Bm, m, kc, c_e, c_h)
    return (Er, Eθ, Ez, Hr, Hθ, Hz)
end

function tm_coax_fields(r::AbstractArray{T, N}, θ::AbstractArray{T, N}, a, b, m, n, f, μᵣ, εᵣ) where {T, N}

    kc = kc_coax(a, b, m, n, :TM)
    β = phase_constant(kc, f, μᵣ, εᵣ)
    Am, Bm = coax_boundary_coeff_tm(m, b, kc)
    c_e, c_h = tm_coefficients(kc, β, f, μᵣ, εᵣ)

    fields = similar(r, NTuple{6, Complex{T}})
    for idx in eachindex(r)
        fields[idx] = tm_coax_fields(r[idx], θ[idx], Am, Bm, m, kc, c_e, c_h)
    end
    return fields
end


"""
    te_coax_fields(r, θ, Am, Bm, m, kc, c_e, c_h)
    te_coax_fields(r, θ, a, b, m, n, f, μᵣ, εᵣ)

# Arguments
- `r`: r polar coordinate.
- `θ`: angle in polar coordinates.
- `Am`: First value returned by `coax_boundary_coeff_te`. Is always 1.0
- `Bm`: Second value retuened by `coax_boundary_coeff_te`.
- `m`:  The azimuthal mode number.
- `n`: radial mode number.
- `kc`: cutoff wavenumber.
- `c_e`: im*ω*μ/kc^2
- `c_h`: -im*β/kc^2
- `a`: outer radius.
- `b`: inner radius
- `f`: frequency.
- `μᵣ`: relative permeability
- `εᵣ`: relative permittivity

# Returns
- `(Er, Eθ, Ez, Hr, Hθ, Hz)`

Calculate the electric and magnetic field components in polar coordinates of a TE_{m,n} mode given the dimensions of the circular section with radius r. 
It is assumed that the center of the circular section is at (0,0).

"""
function te_coax_fields(r, θ, Am, Bm, m, kc, c_e, c_h)

    ∂ψᵢ, ∂ψⱼ, ψₖ = coax_modal_f(r, θ, kc, Am, Bm, m)

    Er = -c_e * ∂ψⱼ
    Eθ = +c_e * ∂ψᵢ
    Ez = zero(Er)
    Hr = -c_h * ∂ψᵢ
    Hθ = -c_h * ∂ψⱼ
    Hz =  -im *  ψₖ

    return (Er, Eθ, Ez, Hr, Hθ, Hz)
end

function te_coax_fields(r, θ, a, b, m, n, f, μᵣ, εᵣ)

    kc = kc_coax(a, b, m, n, :TE)
    β = phase_constant(kc, f, μᵣ, εᵣ)
    Am, Bm = coax_boundary_coeff_te(m, b, kc)
    c_e, c_h = te_coefficients(kc, β, f, μᵣ, εᵣ)

    (Er, Eθ, Ez, Hr, Hθ, Hz) = te_coax_fields(r, θ, Am, Bm, m, kc, c_e, c_h)
    return (Er, Eθ, Ez, Hr, Hθ, Hz)
end

function te_coax_fields(r::AbstractArray{T, N}, θ::AbstractArray{T, N}, a, b, m, n, f, μᵣ, εᵣ) where {T, N}

    kc = kc_coax(a, b, m, n, :TE)
    β = phase_constant(kc, f, μᵣ, εᵣ)
    Am, Bm = coax_boundary_coeff_te(m, b, kc)
    c_e, c_h = te_coefficients(kc, β, f, μᵣ, εᵣ)

    fields = similar(r, NTuple{6, Complex{T}})
    for idx in eachindex(r)
        fields[idx] = te_coax_fields(r[idx], θ[idx], Am, Bm, m, kc, c_e, c_h)
    end
    return fields
end



"""
    tem_coax_fields(r, a, b, μᵣ, εᵣ)

Calculate the electric and magnetic field components of the TEM mode in a coaxial waveguide.

# Arguments
- `r`: Radial coordinate
- `a`: Outer radius
- `b`: Inner radius
- `μᵣ`: Relative permeability
- `εᵣ`: Relative permittivity

# Returns
- `(Er, Eθ, Ez, Hr, Hθ, Hz)`
"""
function tem_coax_fields(r, a, b, μᵣ, εᵣ)

    ε = εᵣ * _εₒ
    μ = μᵣ * _μₒ
    η = sqrt(μ / ε)
    Z0 = sqrt(μ / ε) * log(a / b) / (2π)

    Er = sqrt(2 * Z0) / (r * log(a / b))
    #Er = - 1 / (sqrt(μ/ε)) * 1/r
    Eθ = zero(Er)
    Ez = zero(Er)
    Hr = zero(Er)
    Hθ = Er / η
    Hz = zero(Er)

    return (Er, Eθ, Ez, Hr, Hθ, Hz)
end

tem_normalization_coax() = 1.0

"""
    tm_normalization_coax(a, b, m, kc, β, f, μᵣ, εᵣ)

Normalization factor for TM modes to achieve unit power.

The expression can be derived by integrating the Poynting vector over the cross-section of the guide.

# Arguments
- `a`: Outer radius
- `b`: Inner radius
- `m`: Azimuthal mode index
- `kc`: Cutoff wavenumber
- `β`: Phase constant
- `f`: Frequency in Hz
- `μᵣ`: Relative permeability
- `εᵣ`: Relative permittivity
"""
function tm_normalization_coax(a, b, m, kc, β, f, μᵣ, εᵣ)
    ω = 2π * f
    ε = εᵣ * _εₒ
    μ = μᵣ * _μₒ
    Am = 1.0
    Bm = -besselj(m, kc * b) / bessely(m, kc * b)

    F(m, r, kc, Am, Bm) = Am * besselj(m, r * kc) + Bm * bessely(m, r * kc)
    F_prime(m, r, kc, Am, Bm) = 1 / kc * (Am * besselj_prime(m, r * kc) + Bm * bessely_prime(m, r * kc))
    parte1 = (a * F(m, a, kc, Am, Bm) * F_prime(m, a, kc, Am, Bm)) - (b * F(m, b, kc, Am, Bm) * F_prime(m, b, kc, Am, Bm))

    # Lommel integral
    xa = kc * a
    xb = kc * b
    parte2 = Am^2 * lommel_integral(besselj, besselj, m, xa)
    parte2 += Bm^2 * lommel_integral(bessely, bessely, m, xa)
    parte2 += 2 * Am * Bm * lommel_integral(besselj, bessely, m, xa)

    parte2 -= Am^2 * lommel_integral(besselj, besselj, m, xb)
    parte2 -= Bm^2 * lommel_integral(bessely, bessely, m, xb)
    parte2 -= 2 * Am * Bm * lommel_integral(besselj, bessely, m, xb)

    Nmn_TM = (ω * ε * β / (kc^4)) * π * (parte1 + 1 * parte2)
    if m == 0
        Nmn_TM = Nmn_TM * 2
    end

    return sqrt(2 * 1 / Nmn_TM)
end


"""
    te_normalization_coax(a, b, m, kc, β, f, μᵣ, εᵣ)

Normalization factor for TE modes to achieve unit power.

The expression can be derived by integrating the Poynting vector over the cross-section of the guide.

# Arguments
- `a`: Outer radius
- `b`: Inner radius
- `m`: Azimuthal mode index
- `kc`: Cutoff wavenumber
- `β`: Phase constant
- `f`: Frequency in Hz
- `μᵣ`: Relative permeability
- `εᵣ`: Relative permittivity
"""
function te_normalization_coax(a, b, m, kc, β, f, μᵣ, εᵣ)
    ω = 2π * f
    ε = εᵣ * _εₒ
    μ = μᵣ * _μₒ
    Am = 1.0
    Bm = -besselj_prime(m, kc * b) / bessely_prime(m, kc * b)

    F(m, r, kc, Am, Bm) = Am * besselj(m, r * kc) + Bm * bessely(m, r * kc)
    F_prime(m, r, kc, Am, Bm) = 1 / kc * (Am * besselj_prime(m, r * kc) + Bm * bessely_prime(m, r * kc))
    parte1 = (a * F(m, a, kc, Am, Bm) * F_prime(m, a, kc, Am, Bm)) - (b * F(m, b, kc, Am, Bm) * F_prime(m, b, kc, Am, Bm))

    # Lommel integral
    xa = kc * a
    xb = kc * b
    parte2 = Am^2 * lommel_integral(besselj, besselj, m, xa)
    parte2 += Bm^2 * lommel_integral(bessely, bessely, m, xa)
    parte2 += 2 * Am * Bm * lommel_integral(besselj, bessely, m, xa)

    parte2 -= Am^2 * lommel_integral(besselj, besselj, m, xb)
    parte2 -= Bm^2 * lommel_integral(bessely, bessely, m, xb)
    parte2 -= 2 * Am * Bm * lommel_integral(besselj, bessely, m, xb)

    Nmn_TE = (ω * μ * β / (kc^4)) * π * (parte1 + 1 * parte2)
    if m == 0
        Nmn_TE = Nmn_TE * 2
    end

    return sqrt(2 * 1 / Nmn_TE)
end
