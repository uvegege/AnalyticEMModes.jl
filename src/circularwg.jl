using Bessels
using LinearAlgebra

# n = 1:5, m = 0:11
const _x¹mn = collect([3.8318 1.8412 3.0542 4.2012 5.3175 6.4155 7.5013 8.5777 9.6474 10.7114 11.7708 12.8264;
    7.0156 5.3315 6.7062 8.0153 9.2824 10.5199 11.7349 12.9324 14.1155 15.2867 16.4479 17.6003;
    10.1735 8.5363 9.9695 11.3459 12.6819 13.9872 15.2682 16.5294 17.7740 19.0046 20.2230 21.4309;
    13.3237 11.7060 13.1704 14.5859 15.9641 17.3129 18.6375 19.9419 21.2291 22.5014 23.7607 25.0085;
    16.4706 14.8636 16.3475 17.7888 19.1960 20.5755 21.9317 23.2681 24.5872 25.8913 27.1820 28.4609]')

const _xmn = collect([2.4049 3.8318 5.1357 6.3802 7.5884 8.7715 9.9361 11.0864 12.2251 13.3543 14.4755 15.5898;
    5.5201 7.0156 8.4173 9.7610 11.0647 12.3386 13.5893 14.8213 16.0378 17.2412 18.4335 19.6160;
    8.6537 10.1735 11.6199 13.0152 14.3726 15.7002 17.0038 18.2876 19.5545 20.8071 22.0470 23.2759;
    11.7915 13.3237 14.7960 16.2235 17.6160 18.9801 20.3208 21.6415 22.9452 24.2339 25.5095 26.7733;
    14.9309 16.4706 17.9598 19.4094 20.8269 22.2178 23.5861 24.9349 26.2668 27.5838 28.8874 30.1791]')

kc_te_cwg(r, m, n) = _x¹mn[m+1, n] / r
kc_tm_cwg(r, m, n) = _xmn[m+1, n] / r


# https://dlmf.nist.gov/10.21#vi
"""
    _jvm(a, μ)

McMahon's Asymptotic Expansions for Large Zero of Bessel function

Reference: # https://dlmf.nist.gov/10.21#vi
"""
_jvm(a, μ) = a - (μ - 1)/(8*a) - (4*μ - 1)*(7*μ - 31)/(3*(8*a)^3) - (32*(μ-1)*(83*μ^2-982*μ+3779))/(15*(8*a)^5) - 64*(μ - 1)*(6949*μ^3 - 153855*μ^2 + 1585743*μ - 6277237)/(105*(8*a)^7)

"""
    _jvm(a, μ)
    
McMahon's Asymptotic Expansions for Large Zero of Bessel function first derivative.

Reference: # https://dlmf.nist.gov/10.21#vi

"""
_j1vm(b, μ) = b - (μ + 3)/(8*b) - 4*(7*μ^2 + 82*μ - 9)/(3*(8*b)^3) - (32*(83*μ^3 + 2075*μ^2 - 3039*μ + 3537))/(15*(8*b)^5) - 64*(6949*μ^4 + 296492*μ^3 - 1248002*μ^2 + 7414380*μ - 5853627)/(105*(8*b)^7)


"""
    asymtotic_zero_tm(m, n)

McMahon's Asymptotic Expansions for Large Zero of Bessel function.

Reference: # https://dlmf.nist.gov/10.21#vi
"""
function asymtotic_zero_tm(m, n)
    a = π * (n + m/2 - 1/4)
    μ = 4*m^2
    return _jvm(a, μ)
end


"""
    asymtotic_zero_te(m, n)
    
McMahon's Asymptotic Expansions for Large Zero of Bessel function first derivative.

Reference: # https://dlmf.nist.gov/10.21#vi
"""
function asymtotic_zero_te(m, n)
    b = π * (n + m/2 - 3/4)
    μ = 4*m^2
    return _j1vm(b, μ)
end


"""
    Δmax_asymptotic(kind, m, n)

An arbitrary criterion that may or may not be useful; it should be tested further and removed if it is not necessary.
"""
function Δmax_asymptotic(kind, m, n)
    if kind == :TM
        s1 = π * (n + m/2 - 1/4)
        s2 = π * ((n+1) + m/2 - 1/4)
    else
        s1 =  π * (n + m/2 - 3/4)
        s2 =  π * ((n+1) + m/2 - 3/4)
    end
    spacing = s2 - s1
    return min(0.6 * spacing, π/2)
end

"""
    kc_cwg(r, m, n, T; tol = 1e-9, max_iters = 10)

Calculates the cutoff wave number of the m,n mode for a circular waveguide with radius r.
Uses tabulated data or asymptotic expansions of Bessel functions to calculate the zero.

# Arguments
- `r`: Radius of the circular waveguide
- `m`: Azimuthal mode index
- `n`: Radial mode index
- `T`: Mode type (`:TE` or `:TM`)
- `tol`: Tolerance for zero-finding (default: 1e-9)
- `max_iters`: Maximum iterations for Newton-Raphson (default: 10)

# Example
```julia
kc_cwg(1.0, 2, 1, :TE)
kc_cwg(1.0, 2, 1, :TM)
```
"""
function kc_cwg(r, m, n, T; tol = 1e-9, max_iters = 10)

    if !(T == :TE || T ==:TM)
        @error "$T not suported on circular waveguides."
    end
    
    zval = zero(r)
    if (0 <= m <= 11) && (1 <= n <= 5) 
        if T == :TE
            zval = kc_te_cwg(r, m, n)*r
        else
            zval = kc_tm_cwg(r, m, n)*r
        end
    else
        if n >= 0 # Asymptotic zeros. It works well when n >= 5. 
            if T == :TE
                zval = asymtotic_zero_te(m, n)
            else
                zval = asymtotic_zero_tm(m, n)
            end
        end
    end
    
    iters = 0
    δ = Δmax_asymptotic(T, m, n)
    if T == :TE
        x = zval
        f = besselj_prime(m, x)
        while (abs(f) > tol) && (iters < max_iters)
            Δ = f / (-(1/x) * f - ((x^2 - m^2) / x^2) * besselj(m, x))
            if abs(Δ) > δ   
                Δ = sign(Δ) * δ
            end
            x -= Δ
            iters += 1
            f = besselj_prime(m, x)
        end
        zval = x
    else
        x = zval
        f = besselj(m, x)
        while (abs(f) > tol) && (iters < max_iters)
            Δ = f / besselj_prime(m, x)
            if abs(Δ) > δ  
                Δ = sign(Δ) * δ
            end
            x -= Δ 
            iters += 1
            f = besselj(m, x)
        end
        zval = x
    end

    return zval/r
end


"""
    m_over_r_Jm(m, k, r; x_th=1e-3, tol=1e-13, maxs=10)

Computes the product m/r * besselj(m, kc*r). Returns (m / r) * besselj(m, x) except if `iszero(r)`.
"""
function m_over_r_Jm(m, k, r)
    x = k * r
    if iszero(x)
        return m == 0 ? 0.0 : (m == 1 ? k / 2 : 0.0)
    else
        return (m / r) * besselj(m, x)
    end
end



"""
    cwg_modal_f(r, θ, kc, m)

Modal functions for a TE/TM mode on a circular waveguide.

# Arguments
- `r`: Radial coordinate
- `θ`: Angular coordinate
- `kc`: Cutoff wavenumber
- `m`: Azimuthal mode index

# Returns
A tuple (∂Fz/∂ρ, ∂Fz/∂θ, Fz) for a given (r, θ).
"""
function cwg_modal_f(r, θ, kc, m)
    smθ, cmθ = sincos(m * θ)
    Jm_kcr = besselj(m, kc*r)
    J¹m_kcr = kc * besselj_prime(m,  kc*r) # d(Jm(kc*r))/dr = = J'm(kc*r) * kc 
    m_over_rJm = m_over_r_Jm(m, kc, r)

    ∂ψᵢ = J¹m_kcr * cmθ # ∂Fz/∂r
    ∂ψⱼ = m_over_rJm * -smθ # ∂Fz/∂θ
    ψₖ = Jm_kcr * cmθ # Fz
    return (∂ψᵢ, ∂ψⱼ, ψₖ)
end



"""
    tm_cwg_fields(r, θ, m, kc, c_e, c_h)
    tm_cwg_fields(r, θ, a, m, n, f, μᵣ, εᵣ)

Computes the electric and magnetic field components of a TM_{m,n} mode in a circular waveguide.

# Arguments
- `r`: Radial coordinate
- `θ`: Angular coordinate
- `m`: Azimuthal mode index
- `kc`: Cutoff wavenumber, or
- `a`: Radius of the waveguide
- `n`: Radial mode index
- `f`: Frequency in Hz
- `μᵣ`: Relative permeability
- `εᵣ`: Relative permittivity
- `c_e`, `c_h`: Field coupling coefficients (from `tm_coefficients`)

# Returns
`(Er, Eθ, Ez, Hr, Hθ, Hz)` in cylindrical coordinates centered at (0,0).
"""
function tm_cwg_fields(r, θ, m, kc, c_e, c_h)

    ∂ψᵢ, ∂ψⱼ, ψₖ = cwg_modal_f(r, θ, kc, m)

    Er = -c_e * ∂ψᵢ
    Eθ = -c_e * ∂ψⱼ 
    Ez =  -im *  ψₖ 
    Hr = +c_h * ∂ψⱼ 
    Hθ = -c_h * ∂ψᵢ 
    Hz = zero(Hr)

    return (Er, Eθ, Ez, Hr, Hθ, Hz)
end

function tm_cwg_fields(r, θ, a, m, n, f, μᵣ, εᵣ)

    kc = kc_cwg(a, m, n, :TM)
    β = phase_constant(kc, f, μᵣ, εᵣ)
    c_e, c_h = tm_coefficients(kc, β, f, μᵣ, εᵣ)

    (Er, Eθ, Ez, Hr, Hθ, Hz) = tm_cwg_fields(r, θ, m, kc, c_e, c_h)
    return (Er, Eθ, Ez, Hr, Hθ, Hz)
end


function tm_cwg_fields(r::AbstractArray{T, N}, θ::AbstractArray{T, N}, a, m, n, f, μᵣ, εᵣ) where {T, N}

    kc = kc_cwg(a, m, n, :TM)
    β = phase_constant(kc, f, μᵣ, εᵣ)
    c_e, c_h = tm_coefficients(kc, β, f, μᵣ, εᵣ)

    fields = similar(r, NTuple{6, Complex{T}})
    for idx in eachindex(r)
        fields[idx] = tm_cwg_fields(r[idx], θ[idx], m, kc, c_e, c_h)
    end
    return fields
end

#TODO: Vertical / Horizontal polarization (-pi/2)
"""
    te_cwg_fields(r, θ, m, kc, c_e, c_h)
    te_cwg_fields(r, θ, a, m, n, f, μᵣ, εᵣ)

Computes the electric and magnetic field components of a TE_{m,n} mode in a circular waveguide.

# Arguments
- `r`: Radial coordinate
- `θ`: Angular coordinate
- `m`: Azimuthal mode index
- `kc`: Cutoff wavenumber, or
- `a`: Radius of the waveguide
- `n`: Radial mode index
- `f`: Frequency in Hz
- `μᵣ`: Relative permeability
- `εᵣ`: Relative permittivity
- `c_e`, `c_h`: Field coupling coefficients (from `te_coefficients`)

# Returns
`(Er, Eθ, Ez, Hr, Hθ, Hz)` in cylindrical coordinates centered at (0,0).
"""
function te_cwg_fields(r, θ, m, kc, c_e, c_h)

    ∂ψᵢ, ∂ψⱼ, ψₖ = cwg_modal_f(r, θ, kc, m)

    Er = -c_e * ∂ψⱼ
    Eθ = +c_e * ∂ψᵢ
    Ez = zero(Er)
    Hr = -c_h * ∂ψᵢ
    Hθ = -c_h * ∂ψⱼ
    Hz = -im  *  ψₖ

    return (Er, Eθ, Ez, Hr, Hθ, Hz)
end

function te_cwg_fields(r, θ, a, m, n, f, μᵣ, εᵣ)

    kc = kc_cwg(a, m, n, :TE)
    β = phase_constant(kc, f, μᵣ, εᵣ)
    c_e, c_h = te_coefficients(kc, β, f, μᵣ, εᵣ)

    (Er, Eθ, Ez, Hr, Hθ, Hz) = te_cwg_fields(r, θ, m, kc, c_e, c_h)
    return (Er, Eθ, Ez, Hr, Hθ, Hz)
end

function te_cwg_fields(r::AbstractArray{T, N}, θ::AbstractArray{T, N}, a, m, n, f, μᵣ, εᵣ) where {T, N}

    kc = kc_cwg(a, m, n, :TE)
    β = phase_constant(kc, f, μᵣ, εᵣ)
    c_e, c_h = te_coefficients(kc, β, f, μᵣ, εᵣ)

    fields = similar(r, NTuple{6, Complex{T}})
    for idx in eachindex(r)
        fields[idx] = te_cwg_fields(r[idx], θ[idx], m, kc, c_e, c_h)
    end
    return fields
end



"""
    te_normalization_cwg(a, m, kc, β, f, μᵣ, εᵣ)

Normalization factor for TE modes to achieve unit power.

The expression can be derived by integrating the Poynting vector over the cross-section of the guide.

# Arguments
- `a`: Radius of the waveguide
- `m`: Azimuthal mode index
- `kc`: Cutoff wavenumber
- `β`: Phase constant
- `f`: Frequency in Hz
- `μᵣ`: Relative permeability
- `εᵣ`: Relative permittivity
"""
function te_normalization_cwg(a, m, kc, β, f, μᵣ, εᵣ)
    ω = 2π * f
    ε = εᵣ * _εₒ
    μ = μᵣ * _μₒ
    χ = kc * a
    Am = 1.0
    parte1 = χ * besselj(m, χ) * 1/kc * (besselj_prime(m, χ))  
    parte2 = Am^2 * lommel_integral(besselj, besselj, m, χ)

    Nmn_TE = (ω * μ * β / (kc^4) ) * π * (parte1 + 1 * parte2)
    if m == 0
        Nmn_TE = Nmn_TE * 2
    end

    return sqrt(2*1/Nmn_TE)
end


"""
    tm_normalization_cwg(a, m, kc, β, f, μᵣ, εᵣ)

Normalization factor for TM modes to achieve unit power.

The expression can be derived by integrating the Poynting vector over the cross-section of the guide.

# Arguments
- `a`: Radius of the waveguide
- `m`: Azimuthal mode index
- `kc`: Cutoff wavenumber
- `β`: Phase constant
- `f`: Frequency in Hz
- `μᵣ`: Relative permeability
- `εᵣ`: Relative permittivity
"""
function tm_normalization_cwg(a, m, kc, β, f, μᵣ, εᵣ)
    ω = 2π * f
    ε = εᵣ * _εₒ
    μ = μᵣ * _μₒ
    χ = kc * a
    Am = 1.0
    parte1 = χ * besselj(m, χ) * 1/kc * (besselj_prime(m, χ))  
    parte2 = Am^2 * lommel_integral(besselj, besselj, m, χ)

    Nmn_TM = (ω * ε * β / (kc^4) ) * π * (parte1 + 1 * parte2)
    if m == 0
        Nmn_TM = Nmn_TM * 2
    end

    return sqrt(2*1/Nmn_TM)
end



"""
    metric_and_unit_cylindrical(r, ϕ)

Computes the **metric coefficients (scale factors)** and the **Cartesian components of the unit vectors** (êᵣ, êᵩ) for the 2D cylindrical (polar) coordinate system (r, ϕ).

This function is used for transforming differential operators (like the gradient or Laplacian) and vectors between the cylindrical and Cartesian systems.

# Arguments
- `r`: The radial coordinate (distance from the z-axis).
- `ϕ`: The azimuthal coordinate (angle with respect to the positive x-axis, in radians).

# Returns
A tuple containing the six components:
- `h_r`: Metric coefficient in the radial direction . Always **1.0**.
- `h_ϕ`: Metric coefficient in the azimuthal direction . Equal to `r`.
- `e_r_x`: X-component of the radial unit vector 
- `e_r_y`: Y-component of the radial unit vector 
- `e_ϕ_x`: X-component of the azimuthal unit vector 
- `e_ϕ_y`: Y-component of the azimuthal unit vector
"""
function metric_and_unit_cylindrical(r, ϕ)

    h_r = 1.0
    h_ϕ = r
    e_r_x = cos(ϕ)
    e_r_y = sin(ϕ)
    e_ϕ_x = -sin(ϕ)
    e_ϕ_y =  cos(ϕ)

    return h_r, h_ϕ, e_r_x, e_r_y, e_ϕ_x, e_ϕ_y
end
