using MathieuF
#using LinearAlgebra
include("./mathieu_functions.jl")

# References:
# Modes of Elliptical Waveguides; A Correction. D.A Goldberg, L.J. Laslett and R.A. Rimmer 
# Wave Propagation in Hollow Conducting Elliptical Waveguides. Jan G. Kretzschmar.

"""
    Ce_m(m, q, z)

Even Radial Mathieu function and its derivative. It's equivalent to the Modified Mathieu Function of first kind of order `m` with parameter `q` evaluated at `z`.

# Arguments
- `m`: Order of the Mathieu function
- `q`: Parameter related to the geometry (q = (kc*ρ)²/4)
- `z`: Coordinate (typically ξ)
"""
Ce_m(m, q, z) = mathieu_Mce(1, m, q, z)
"""
    Se_m(m, q, z)

Odd Radial Mathieu function and its derivative. It's equivalent to the Modified Mathieu Function of first kind of order `m` with parameter `q` evaluated at `z`.

# Arguments
- `m`: Order of the Mathieu function
- `q`: Parameter related to the geometry (q = (kc*ρ)²/4)
- `z`: Coordinate (typically ξ)
"""
Se_m(m, q, z) = mathieu_Mse(1, m, q, z)
"""
    ce_m(m, q, z)

Even angular Mathieu function of order `m` with parameter `q` evaluated at `z` and its derivative.

# Arguments
- `m`: Order of the Mathieu function
- `q`: Parameter related to the geometry (q = (kc*ρ)²/4)
- `z`: Coordinate (typically η)
"""
ce_m(m::Int, q, z) = mathieu_ce(m, q, z)
"""
    se_m(m, q, z)

Odd angular Mathieu function of order `m` with parameter `q` evaluated at `z` and its derivative.

# Arguments
- `m`: Order of the Mathieu function
- `q`: Parameter related to the geometry (q = (kc*ρ)²/4)
- `z`: Coordinate (typically η)
"""
se_m(m::Int, q, z) = mathieu_se(m, q, z)


"""
    elliptic_modal_f(ξ, η, m, q, coeff)

Modal functions for a TE/TM mode on a elliptical waveguide. Return (∂Fz/∂ξ, ∂Fz/∂η, Fz)
"""
function elliptic_modal_f(ξ, η, m, q, coeff, even)

    U, U¹ = even ? Mce_kernel(1, m, coeff, q, ξ) : Mse_kernel(1, m, coeff, q, ξ)
    V, V¹ = even ? ce_kernel(m, coeff, q, η) : se_kernel(m, coeff, q, η)

    ∂ψᵢ = U¹ * V 
    ∂ψⱼ = U  * V¹ 
    ψₖ =  U  * V

    return (∂ψᵢ, ∂ψⱼ, ψₖ)
end

"""
    te_ewg_fields(ξ, η, ρ, m, even, coeff, q, c_e, c_h)
    te_ewg_fields(ξ, η, a, b, m, n, even, f, μᵣ, εᵣ)

Computes the electric and magnetic field components of a TE_{c/s,m,n} mode in an elliptical waveguide.

# Arguments
- `ξ`: Radial elliptical coordinate
- `η`: Angular elliptical coordinate
- `ρ`: Focal distance √(a²-b²), or
- `a`: Semi-major axis
- `b`: Semi-minor axis
- `m`: Mode order
- `n`: Radial mode index (used with a, b)
- `even`: `true` for even modes (c), `false` for odd modes (s)
- `coeff`: Mathieu function coefficients
- `q`: Parameter q = (kc*ρ)²/4
- `f`: Frequency in Hz
- `μᵣ`: Relative permeability
- `εᵣ`: Relative permittivity
- `c_e`, `c_h`: Field coupling coefficients

# Returns
`(Eξ, Eη, Ez, Hξ, Hη, Hz)` in elliptical coordinates centered at (0,0).
"""
function te_ewg_fields(ξ, η, ρ, m, even, coeff, q, c_e, c_h)

    ∂ψᵢ, ∂ψⱼ, ψₖ = elliptic_modal_f(ξ, η, m, q, coeff, even)
    h = ρ * sqrt(sinh(ξ)^2 + sin(η)^2)

    Eξ = -c_e/h * ∂ψⱼ
    Eη = +c_e/h * ∂ψᵢ
    Ez = zero(Eξ)
    Hξ = -c_h/h * ∂ψᵢ
    Hη = -c_h/h * ∂ψⱼ
    Hz =   -im  *  ψₖ

    return (Eξ, Eη, Ez, Hξ, Hη, Hz)
end

function te_ewg_fields(ξ, η, a, b, m, n, even, f, μᵣ, εᵣ)
    ρ = sqrt(a^2-b^2)
    kc = kc_ewg(a, b, m, n, even ,:TE)
    q = (kc*ρ)^2/4
    c = even ? MathieuCharA(m, q) : MathieuCharB(m, q)
    coeff = even ? mathieu_a_coeff(m, q, c, 100) : mathieu_b_coeff(m, q, c, 100)
    β = phase_constant(kc, f, μᵣ, εᵣ)
    c_e, c_h = te_coefficients(kc, β, f, μᵣ, εᵣ)
    (Eξ, Eη, Ez, Hξ, Hη, Hz) = te_ewg_fields(ξ, η, ρ, m, even, coeff, q, c_e, c_h)
    return (Eξ, Eη, Ez, Hξ, Hη, Hz)
end

function te_ewg_fields(ξ::AbstractArray{T, N}, η::AbstractArray{T, N}, a, b, m, n, even, f , μᵣ, εᵣ) where {T, N}
    ρ = sqrt(a^2-b^2)
    kc = kc_ewg(a, b, m, n, even ,:TE)
    q = (kc*ρ)^2/4
    c = even ? MathieuCharA(m, q) : MathieuCharB(m, q)
    coeff = even ? mathieu_a_coeff(m, q, c, 100) : mathieu_b_coeff(m, q, c, 100)
    β = phase_constant(kc, f, μᵣ, εᵣ)
    c_e, c_h = te_coefficients(kc, β, f, μᵣ, εᵣ)
    fields = similar(ξ, NTuple{6, Complex{T}})
    for idx in eachindex(ξ)
        fields[idx] = te_ewg_fields(ξ[idx], η[idx], ρ, m, even, coeff, q, c_e, c_h)
    end
    return fields
end


"""
    tm_ewg_fields(ξ, η, ρ, m, even, coeff, q, c_e, c_h)
    tm_ewg_fields(ξ, η, a, b, m, n, even, f, μᵣ, εᵣ)

Computes the electric and magnetic field components of a TM_{c/s,m,n} mode in an elliptical waveguide.

# Arguments
- `ξ`: Radial elliptical coordinate
- `η`: Angular elliptical coordinate
- `ρ`: Focal distance √(a²-b²), or
- `a`: Semi-major axis
- `b`: Semi-minor axis
- `m`: Mode order
- `n`: Radial mode index (used with a, b)
- `even`: `true` for even modes (c), `false` for odd modes (s)
- `coeff`: Mathieu function coefficients
- `q`: Parameter q = (kc*ρ)²/4
- `f`: Frequency in Hz
- `μᵣ`: Relative permeability
- `εᵣ`: Relative permittivity
- `c_e`, `c_h`: Field coupling coefficients

# Returns
`(Eξ, Eη, Ez, Hξ, Hη, Hz)` in elliptical coordinates centered at (0,0).
"""
function tm_ewg_fields(ξ, η, ρ, m, even, coeff, q, c_e, c_h)

    ∂ψᵢ, ∂ψⱼ, ψₖ = elliptic_modal_f(ξ, η, m, q, coeff, even)
    h = ρ * sqrt(sinh(ξ)^2 + sin(η)^2)

    Eξ = -c_e/h * ∂ψᵢ
    Eη = -c_e/h * ∂ψⱼ
    Ez =   -im  *  ψₖ 
    Hξ = +c_h/h * ∂ψⱼ
    Hη = -c_h/h * ∂ψᵢ
    Hz = zero(Hξ)

    return (Eξ, Eη, Ez, Hξ, Hη, Hz)
end

function tm_ewg_fields(ξ, η, a, b, m, n, even, f, μᵣ, εᵣ)
    ρ = sqrt(a^2-b^2)
    kc = kc_ewg(a, b, m, n, even ,:TM)
    q = (kc*ρ)^2/4
    c = even ? MathieuCharA(m, q) : MathieuCharB(m, q)
    coeff = even ? mathieu_a_coeff(m, q, c, 100) : mathieu_b_coeff(m, q, c, 100)
    β = phase_constant(kc, f, μᵣ, εᵣ)
    c_e, c_h = tm_coefficients(kc, β, f, μᵣ, εᵣ)
    (Eξ, Eη, Ez, Hξ, Hη, Hz) = tm_ewg_fields(ξ, η, ρ, m, even, coeff, q, c_e, c_h)
    return (Eξ, Eη, Ez, Hξ, Hη, Hz)
end

function tm_ewg_fields(ξ::AbstractArray{T, N}, η::AbstractArray{T, N}, a, b, m, n, even, f , μᵣ, εᵣ) where {T, N}
    ρ = sqrt(a^2-b^2)
    kc = kc_ewg(a, b, m, n, even ,:TM)
    q = (kc*ρ)^2/4
    c = even ? MathieuCharA(m, q) : MathieuCharB(m, q)
    coeff = even ? mathieu_a_coeff(m, q, c, 100) : mathieu_b_coeff(m, q, c, 100)
    β = phase_constant(kc, f, μᵣ, εᵣ)
    c_e, c_h = tm_coefficients(kc, β, f, μᵣ, εᵣ)

    fields = similar(ξ, NTuple{6, Complex{T}})
    for idx in eachindex(ξ)
        fields[idx] = tm_ewg_fields(ξ[idx], η[idx], ρ, m, even, coeff, q, c_e, c_h)
    end

    return fields
end


#TODO: Check m >> 1
"""
    kc_ewg(a, b, m, n, even, T)

Computes the cutoff wavenumber `k_c` for the (m, n)-th mode in an **elliptical waveguide** (EWG).

This function determines the characteristic values of the **Mathieu functions** corresponding to
the given mode symmetry (`even`) and type (`T = :TE` or `:TM`).  
The cutoff wavenumber is derived from the eigenvalues `qₘₙ` obtained by locating the zeros of the
appropriate Mathieu function or its derivative.

# Arguments
- `a`, `b`: Semi-major and semi-minor axes of the ellipse.  
- `m`, `n`: Mode indices.  
- `even::Bool`: Whether the mode is even (`true`) or odd (`false`).  
- `T::Symbol`: Mode type. Either `:TE` or `:TM`.

# Returns
- `kc::Float64`: Cutoff wavenumber for the specified mode.

# Notes
- For TM modes, the zeros of the Mathieu functions Ceₘ or Seₘ are used.
- For TE modes, the zeros of their derivatives Ceₘ′ or Seₘ′ are used.
"""
function kc_ewg(a, b, m, n, even, T)
    ρ = sqrt(a^2 - b^2)
    e = ρ / a
    ξ = acosh(1 / e)
    if T == :TM
        if even
            return sqrt(4*find_n_zero_ewg(q->Ce_m(m, q, ξ), m, n, T))/ρ
        else
            return sqrt(4*find_n_zero_ewg(q->Se_m(m, q, ξ), m, n, T))/ρ
        end
    elseif T == :TE
        if even
            return sqrt(4*find_n_zero_ewg(q->Ce_m(m, q, ξ)[2], m, n, T))/ρ
        else
            return sqrt(4*find_n_zero_ewg(q->Se_m(m, q, ξ)[2], m, n, T))/ρ
        end
    else 
        @error "Not supported mode"
    end
end

"""
    metric_and_unit_elliptic(ρ, ξ, η)

Computes the **metric factor** and **unit vectors** for the 2D elliptic coordinate system.
"""
function metric_and_unit_elliptic(ρ, ξ, η)
    s1 = sinh(ξ); s2 = sin(η)
    denom = sqrt(s1^2 + s2^2)
    e_xi_x = s1 * cos(η) / denom
    e_xi_y = cosh(ξ) * sin(η) / denom
    e_eta_x = -cosh(ξ) * sin(η) / denom
    e_eta_y =  sinh(ξ) * cos(η) / denom
    h = ρ * denom
    return h, e_xi_x, e_xi_y, e_eta_x, e_eta_y
end


"""
    cart2elliptic(x, y, a, b)

Transform from cartesian coordinates to elliptic coordinates.

# Arguments
- `x`, `y`: Cartesian coordinates
- `a`: Semi-major axis of the ellipse
- `b`: Semi-minor axis of the ellipse

# Returns
`(ξ, η)` elliptic coordinates
"""
function cart2elliptic(x, y, a, b)
    c = sqrt(a^2-b^2)
    B = x^2+y^2-c^2
    p = (-B + sqrt(B^2+4*c^2*y^2))/(2*c^2)
    q = (-B - sqrt(B^2+4*c^2*y^2))/(2*c^2)
    η₀ = asin(sqrt(p))
    if x >= 0 && y >= 0
        η = η₀
    elseif x < 0 && y >= 0
        η = π - η₀
    elseif x <= 0 && y < 0
        η = π + η₀
    else
        η = 2π - η₀
    end
    ξ = 1/2 * log(1 - 2*q + 2*sqrt(q^2-q))
    return ξ, η
end
