"""
    kc_wedge(βz, f, μᵣ, εᵣ)

Computes the radial propagation constant for a wedge waveguide (same as `kc_radial`).

# Arguments
- `βz`: Axial phase constant
- `f`: Frequency in Hz
- `μᵣ`: Relative permeability
- `εᵣ`: Relative permittivity
"""
kc_wedge(βz, f, μᵣ, εᵣ) = kc_radial(βz, f, μᵣ, εᵣ)

"""
    te_wedge_fields(r, ϕ, z, ϕ0, p, Amn, Bmn, kc, β, c_e, c_h)
    te_wedge_fields(r, ϕ, z, h, ϕ0, p, n, Amn, Bmn, f, μᵣ, εᵣ)

Compute the electric and magnetic field components of a TE_{p,n} mode in a wedge waveguide.

# Arguments
- `r`: Radial coordinate
- `ϕ`: Angular coordinate where ϕ ∈ [0, ϕ0]
- `z`: Axial coordinate where z ∈ [0, h]
- `ϕ0`: Wedge angle
- `p`: Angular mode index
- `Amn`, `Bmn`: Wave amplitudes
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
"""
te_wedge_fields(r, ϕ, z, ϕ0, p, Amn, Bmn, kc, β, c_e, c_h)  = te_radial_fields(r, ϕ, z, p*π/ϕ0, Amn, Bmn, kc, β, c_e, c_h)
te_wedge_fields(r, ϕ, z, h, ϕ0, p, n, Amn, Bmn, f, μᵣ, εᵣ) = te_radial_fields(r, ϕ, z, h, p*π/ϕ0, n, Amn, Bmn, f, μᵣ, εᵣ)

"""
    tm_wedge_fields(r, ϕ, z, ϕ0, p, Amn, Bmn, kc, β, c_e, c_h)
    tm_wedge_fields(r, ϕ, z, h, ϕ0, p, n, Amn, Bmn, f, μᵣ, εᵣ)

Compute the electric and magnetic field components of a TM_{p,n} mode in a wedge waveguide.

# Arguments
- `r`: Radial coordinate
- `ϕ`: Angular coordinate where ϕ ∈ [0, ϕ0]
- `z`: Axial coordinate where z ∈ [0, h]
- `ϕ0`: Wedge angle
- `p`: Angular mode index
- `Amn`, `Bmn`: Wave amplitudes
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
"""
tm_wedge_fields(r, ϕ, z, ϕ0, p, Amn, Bmn, kc, β, c_e, c_h)  = tm_radial_fields(r, ϕ, z, p*π/ϕ0, Amn, Bmn, kc, β, c_e, c_h)
tm_wedge_fields(r, ϕ, z, h, ϕ0, p, n, Amn, Bmn, f, μᵣ, εᵣ) = tm_radial_fields(r, ϕ, z, h, p*π/ϕ0, n, Amn, Bmn, f, μᵣ, εᵣ)


"""
    te_normalization_wedge(r, h, ϕ0, p, Amn, Bmn, kc, f, μᵣ, εᵣ)

Normalization factor for TE modes to achieve unit power.

The expression can be derived by integrating the Poynting vector over the cross-section of the guide.

# Arguments
- `r`: Radial coordinate
- `h`: Height of the waveguide
- `ϕ0`: Wedge angle
- `p`: Angular mode index
- `Amn`, `Bmn`: Wave amplitudes
- `kc`: Radial propagation constant
- `f`: Frequency in Hz
- `μᵣ`: Relative permeability
- `εᵣ`: Relative permittivity
"""
function te_normalization_wedge(r, h, ϕ0, p, Amn, Bmn, kc, f, μᵣ, εᵣ)
    ω = 2π * f
    μ = μᵣ * _μₒ
    m = p*π/ϕ0
    F = (Amn*hankelh2(m, kc*r) + Bmn*hankelh1(m, kc * r))
    F_prime = kc * (Amn * hankelh2_prime(m, kc*r)+ Bmn * hankelh1_prime(m, kc*r)) 
    Nmn = 1/2 * im*μ*ω/kc^2 * h/2 * ifelse(m == 0, ϕ0, ϕ0/2) * r * real(im * conj(F) * F_prime)
    return sqrt(1/abs(Nmn))
end


"""
    tm_normalization_wedge(r, h, ϕ0, p, Amn, Bmn, kc, f, μᵣ, εᵣ)

Normalization factor for TM modes to achieve unit power.

The expression can be derived by integrating the Poynting vector over the cross-section of the guide.

# Arguments
- `r`: Radial coordinate
- `h`: Height of the waveguide
- `ϕ0`: Wedge angle
- `p`: Angular mode index
- `Amn`, `Bmn`: Wave amplitudes
- `kc`: Radial propagation constant
- `f`: Frequency in Hz
- `μᵣ`: Relative permeability
- `εᵣ`: Relative permittivity
"""
function tm_normalization_wedge(r, h, ϕ0, p, Amn, Bmn, kc, f, μᵣ, εᵣ)
    ω = 2π * f
    ε = εᵣ * _εₒ
    m = p*π/ϕ0
    F = (Amn*hankelh2(m, kc*r) + Bmn*hankelh1(m, kc * r))
    F_prime = kc * (Amn * hankelh2_prime(m, kc*r)+ Bmn * hankelh1_prime(m, kc*r)) 
    Nmn = 1/2*ε*ω/kc^2 * h/2 * ifelse(m == 0, ϕ0, ϕ0/2) * r * real(im * F * conj(F_prime))
    return sqrt(1/abs(Nmn))
end
