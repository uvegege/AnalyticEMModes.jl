"""
    phase_constant_radial(h, n)

Phase constant (βᶻ) of the **axial component** of the wavenumber in a radial waveguide.

This constant arises from the **PEC boundary conditions** imposed at `z = 0` and `z = h`,
which constrain the field variation along the axial direction as:

Eρ(0 <= ρ <= ∞, 0 <= ϕ <= 2π, z = 0) = Eρ(0 <= ρ <= ∞, 0 <= ϕ <= 2π, z = h) = 0
Eϕ(0 <= ρ <= ∞, 0 <= ϕ <= 2π, z = 0) = Eϕ(0 <= ρ <= ∞, 0 <= ϕ <= 2π, z = h) = 0

# Arguments
- `h`: Height of the radial waveguide
- `n`: Axial mode index

# Reference
Advanced engineering electromagnetics by Constatine A. Balanis. Chapter 9.5.1
"""
phase_constant_radial(h, n) = n*π/h 

"""
    kc_radial(βz, f, μᵣ, εᵣ)

Computes the **radial propagation constant** `kᵨ` (analogous to `kc` in other geometries)
for a **radial waveguide**.

# Arguments
- `βz`: Axial phase constant
- `f`: Frequency in Hz
- `μᵣ`: Relative permeability
- `εᵣ`: Relative permittivity

# Reference
Advanced engineering electromagnetics by Constatine A. Balanis. Chapter 9.5.1
"""
function kc_radial(βz, f, μᵣ, εᵣ)
    μ = μᵣ * _μₒ
    ε = εᵣ * _εₒ
    ω = 2π * f
    k = ω * sqrt(μ * ε)
    v = sqrt(Complex(k^2 - βz^2))
    #if imag(v) != 0
    #    @error "Not implemented for frequencies below cutoff" 
    #end
    return real(v)
end

"""
    cutoff_frequency_radial(βz, μᵣ, εᵣ)

Like `cutoff_frequency(kc, μᵣ, εᵣ)` but using `βz` instead of `kc` for the particular case of the radial waveguide.

# Arguments
- `βz`: Axial phase constant
- `μᵣ`: Relative permeability
- `εᵣ`: Relative permittivity
"""
cutoff_frequency_radial(βz, μᵣ, εᵣ) = cutoff_frequency(βz, μᵣ, εᵣ)


"""
    radial_modal_f(r, ϕ, kc, m, Amn, Bmn)

Modal functions for a TE/TM mode on a radial waveguide.

# Arguments
- `r`: Radial coordinate
- `ϕ`: Angular coordinate
- `kc`: Radial propagation constant
- `m`: Azimuthal mode index
- `Amn`, `Bmn`: Wave amplitudes

# Notes
Outgoing waves: `F⁺z(ρ, ϕ, z) = Amn * H⁽²⁾ₘ(kc * ρ) * cos(mϕ) * sin(nπ/h * z)` → Amn = 1, Bmn = 0

Incoming waves: `F⁻z(ρ, ϕ, z) = Bmn * H⁽¹⁾ₘ(kc * ρ) * cos(mϕ) * sin(nπ/h * z)` → Amn = 0, Bmn = 1
"""
function radial_modal_f(r, ϕ, kc, m, Amn, Bmn)

    smϕ, cmϕ =  sincos(m*ϕ)

    ψ =(Amn*hankelh2(m,kc*r) + Bmn*hankelh1(m, kc*r))
    ψ¹ = kc * (Amn*hankelh2_prime(m,kc*r)+ Bmn*hankelh1_prime(m, kc*r))
    
    ∂ψᵢ = ψ¹ * cmϕ # ∂Fz/dρ
    ∂ψⱼ = m/r * ψ * -smϕ # ∂Fz/dϕ
    ψₖ = ψ * cmϕ

    return (∂ψᵢ, ∂ψⱼ, ψₖ)
end


"""
    te_radial_fields(r, ϕ, z, m, Amn, Bmn, kc, β, c_e, c_h)
    te_radial_fields(r, ϕ, z, h, m, n, Amn, Bmn, f, μᵣ, εᵣ)

Compute the electric and magnetic field components of a TE_{m,n} mode in a radial waveguide.

# Arguments
- `r`: Radial coordinate
- `ϕ`: Angular coordinate
- `z`: Axial coordinate where z ∈ [0, h]
- `m`: Azimuthal mode index
- `Amn`, `Bmn`: Wave amplitudes (outgoing/incoming)
- `kc`: Radial propagation constant, or
- `h`: Height of the waveguide
- `n`: Axial mode index
- `f`: Frequency in Hz
- `μᵣ`: Relative permeability
- `εᵣ`: Relative permittivity
- `β`: Axial phase constant
- `c_e`, `c_h`: Field coupling coefficients

# Returns
`(Er, Eϕ, Ez, Hr, Hϕ, Hz)`

# Notes
Outgoing: Amn = 1, Bmn = 0; Incoming: Amn = 0, Bmn = 1
"""
function te_radial_fields(r, ϕ, z, m, Amn, Bmn, kc, β, c_e, c_h)

    ∂ψᵢ, ∂ψⱼ, ψₖ = radial_modal_f(r, ϕ, kc, m, Amn, Bmn)
    snz, cnz = sincos(β*z)

    Er = -c_e * ∂ψⱼ * snz
    Eϕ = +c_e * ∂ψᵢ * snz
    Ez = zero(Er)
    Hr = -im * c_h * ∂ψᵢ * cnz
    Hϕ = -im * c_h * ∂ψⱼ * cnz
    Hz =  -im *  ψₖ * snz

    return (Er, Eϕ, Ez, Hr, Hϕ, Hz)
end

function te_radial_fields(r, ϕ, z, h, m, n, Amn, Bmn, f, μᵣ, εᵣ)
    β = phase_constant_radial(h, n)
    kc = kc_radial(β, f, μᵣ, εᵣ)
    c_e, c_h = te_coefficients(kc, β, f, μᵣ, εᵣ)

    (Er, Eϕ, Ez, Hr, Hϕ, Hz) = te_radial_fields(r, ϕ, z, m, Amn, Bmn, kc, β, c_e, c_h)
    return (Er, Eϕ, Ez, Hr, Hϕ, Hz)
end

function te_radial_fields(r::AbstractArray{T, N}, ϕ::AbstractArray{T, N}, z::AbstractArray{T, N}, h, m, n, Amn, Bmn, f, μᵣ, εᵣ) where {T, N}
    β = phase_constant_radial(h, n)
    kc = kc_radial(β, f, μᵣ, εᵣ)
    c_e, c_h = te_coefficients(kc, β, f, μᵣ, εᵣ)

    fields = similar(r, NTuple{6, Complex{T}})
    for idx in eachindex(r)
        fields[idx] = te_radial_fields(r[idx], ϕ[idx], z[idx], m, Amn, Bmn, kc, β, c_e, c_h)
    end
    return fields
end


"""
    te_zw(m, r, kc, Amn, Bmn, f, μᵣ)

Impedance of the wave in the ρ direction for a TE mode.

# Arguments
- `m`: Azimuthal mode index
- `r`: Radial coordinate
- `kc`: Radial propagation constant
- `Amn`, `Bmn`: Wave amplitudes
- `f`: Frequency in Hz
- `μᵣ`: Relative permeability
"""
function te_zw(m, r, kc, Amn, Bmn, f, μᵣ)
    μ = μᵣ * _μₒ
    ω = 2*π*f
    return im*ω*μ/kc * (Amn*hankelh2_prime(m,kc*r)+ Bmn*hankelh1_prime(m, kc*r))/(Amn*hankelh2(m, kc*r)+ Bmn*hankelh1(m, kc*r))
end




"""
    tm_radial_fields(r, ϕ, z, m, Amn, Bmn, kc, β, c_e, c_h)
    tm_radial_fields(r, ϕ, z, h, m, n, Amn, Bmn, f,  μᵣ, εᵣ )

Compute the electric and magnetic field components of a TM_{m,n} mode given the dimensions radial waveguide.
It is assumed that the radial section is centered at (0, 0, z) with height z ∈ [0, h].

"""
function tm_radial_fields(r, ϕ, z, m, Amn, Bmn, kc, β, c_e, c_h)

    ∂ψᵢ, ∂ψⱼ, ψₖ = radial_modal_f(r, ϕ, kc, m, Amn, Bmn)
    snz, cnz = sincos(β*z)

    Er = +im * c_e * ∂ψᵢ * snz
    Eϕ = +im * c_e * ∂ψⱼ * snz
    Ez =  -im *  ψₖ * cnz
    Hr = +c_h * ∂ψⱼ * cnz
    Hϕ = -c_h * ∂ψᵢ * cnz
    Hz = zero(Er)

    return (Er, Eϕ, Ez, Hr, Hϕ, Hz)
end

function tm_radial_fields(r, ϕ, z, h, m, n, Amn, Bmn, f, μᵣ, εᵣ)
    β = phase_constant_radial(h, n)
    kc = kc_radial(β, f, μᵣ, εᵣ)
    c_e, c_h = tm_coefficients(kc, β, f, μᵣ, εᵣ)

    (Er, Eϕ, Ez, Hr, Hϕ, Hz) = tm_radial_fields(r, ϕ, z, m, Amn, Bmn, kc, β, c_e, c_h)
    return (Er, Eϕ, Ez, Hr, Hϕ, Hz)
end

function tm_radial_fields(r::AbstractArray{T, N}, ϕ::AbstractArray{T, N}, z::AbstractArray{T, N}, h, m, n, Amn, Bmn, f, μᵣ, εᵣ) where {T, N}
    β = phase_constant_radial(h, n)
    kc = kc_radial(β, f, μᵣ, εᵣ)
    c_e, c_h = tm_coefficients(kc, β, f, μᵣ, εᵣ)

    fields = similar(r, NTuple{6, Complex{T}})
    for idx in eachindex(r)
        fields[idx] = tm_radial_fields(r[idx], ϕ[idx], z[idx], m, Amn, Bmn, kc, β, c_e, c_h)
    end
    return fields
end

"""
    tm_zw(m, r, kc, Amn, Bmn, f, μᵣ)

Impedance of the wave in the ρ direction for a TM mode.
"""
function tm_zw(m, r, kc, Amn, Bmn, f, εᵣ)
    ε = εᵣ * _εₒ
    ω = 2*π*f
    return -im*kc/(ω*ε) * (Amn*hankelh2(m,kc*r) + Bmn*hankelh1(m, kc*r))/(Amn*hankelh2_prime(m,kc*r)+ Bmn*hankelh1_prime(m, kc*r))
end

"""
    te_normalization_radial(r, h, m, Amn, Bmn, kc, f, μᵣ, εᵣ)

Normalization factor for TE modes to achieve unit power.

The expression can be derived by integrating the Poynting vector over the cross-section of the guide.

# Arguments
- `r`: Radial coordinate
- `h`: Height of the waveguide
- `m`: Azimuthal mode index
- `Amn`, `Bmn`: Wave amplitudes
- `kc`: Radial propagation constant
- `f`: Frequency in Hz
- `μᵣ`: Relative permeability
- `εᵣ`: Relative permittivity
"""
function te_normalization_radial(r, h, m, Amn, Bmn, kc, f, μᵣ, εᵣ)
    ω = 2π * f
    μ = μᵣ * _μₒ
    F = (Amn*hankelh2(m, kc*r) + Bmn*hankelh1(m, kc * r))
    F_prime = kc * (Amn * hankelh2_prime(m, kc*r)+ Bmn * hankelh1_prime(m, kc*r)) 
    Nmn = 1/2 * im*μ*ω/kc^2 * h/2 * ifelse(m == 0, 2*π, π) * r * real(im * conj(F) * F_prime)
    return sqrt(1/abs(Nmn))
end


"""
    tm_normalization_radial(r, h, m, Amn, Bmn, kc, f, μᵣ, εᵣ)

Normalization factor for TM modes to achieve unit power.

The expression can be derived by integrating the Poynting vector over the cross-section of the guide.

# Arguments
- `r`: Radial coordinate
- `h`: Height of the waveguide
- `m`: Azimuthal mode index
- `Amn`, `Bmn`: Wave amplitudes
- `kc`: Radial propagation constant
- `f`: Frequency in Hz
- `μᵣ`: Relative permeability
- `εᵣ`: Relative permittivity
"""
function tm_normalization_radial(r, h, m, Amn, Bmn, kc, f, μᵣ, εᵣ)
    ω = 2π * f
    ε = εᵣ * _εₒ
    F = (Amn*hankelh2(m, kc*r) + Bmn*hankelh1(m, kc * r))
    F_prime = kc * (Amn * hankelh2_prime(m, kc*r)+ Bmn * hankelh1_prime(m, kc*r)) 
    Nmn = 1/2*ε*ω/kc^2 * h/2 * ifelse(m == 0, 2*π, π) * r * real(im * F * conj(F_prime))
    return sqrt(1/abs(Nmn))
end
