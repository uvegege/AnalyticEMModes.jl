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



"""
    Av_Avp_from_kernels(m, q, coeff; even=true, Nη=4096)

Compute the angular normalization integrals:

    A_v  = ∫₀²π |V(η)|²    dη
    A_vp = ∫₀²π |∂V/∂η|²  dη

where `V(η)` is `ce_m(η; q)` (`even=true`) or `se_m(η; q)` (`even=false`), evaluated
using the precomputed Fourier coefficients `coeff`.

The integrals are computed via the **spectral (periodic) trapezoidal rule** with `Nη`
uniformly-spaced points on [0, 2π). Because the angular Mathieu functions are
C∞-smooth and 2π-periodic, the trapezoidal rule converges exponentially fast and is
effectively equivalent to applying Parseval's theorem to the Fourier-series
representation of the Mathieu functions. For `Nη ≥ 4096` the result is
indistinguishable from machine-precision exact.

# Arguments
- `m`: Azimuthal mode index
- `q`: Mathieu parameter q = (kc·ρ)²/4
- `coeff`: Fourier coefficients from `mathieu_a_coeff` (even) or `mathieu_b_coeff` (odd)
- `even`: `true` for `ce_m` (even), `false` for `se_m` (odd)
- `Nη`: Number of quadrature points (default 4096)

# Returns
`(A_v, A_vp)` as a plain tuple.
"""
function Av_Avp_from_kernels(m::Int, q, coeff; even::Bool=true, Nη::Int=4096)
    m = abs(m)

    # se_0 ≡ 0
    if !even && m == 0
        return (0.0, 0.0)
    end

    kernel = even ? ce_kernel : se_kernel
    Δη = 2π / Nη

    Av = 0.0
    Avp = 0.0

    # trapecio periódico: suma uniforme sin endpoints duplicados
    for j in 0:Nη-1
        η = j * Δη
        V, Vη = kernel(m, coeff, q, η)
        Av  += abs2(V)
        Avp += abs2(Vη)
    end

    Av  *= Δη
    Avp *= Δη
    return Av, Avp
end


elliptic_xi0(a,b) = acosh(a/sqrt(a^2 - b^2))

"""
    elliptic_modal_power_1d(a, b, m, q, coeff, kc, β, f, μᵣ, εᵣ;
                                       pol=:TE, even=true,
                                       Nη=4096, rtol=1e-10, atol=0.0, maxdepth=25)

Time-averaged modal power P (W) for a mode with unit longitudinal amplitude, computed
by integrating the Poynting vector over the elliptical waveguide cross-section.

In elliptic coordinates (ξ, η), the cross-sectional Poynting integral factorizes into
a product of independent 1D integrals:

    P = (1/2) · (ω·σ·β / kc⁴) · (A_v · I_U′ + A_vp · I_U)

where σ = μ for TE modes and σ = ε for TM modes, and

    A_v  = ∫₀²π |V(η)|²    dη ,   A_vp = ∫₀²π |∂V/∂η|²  dη
    I_U  = ∫₀^{ξ₀} |U(ξ)|² dξ ,   I_U′ = ∫₀^{ξ₀} |∂U/∂ξ|² dξ

**Angular integrals (quasi-analytical):** `A_v` and `A_vp` are evaluated via the
spectral trapezoidal rule applied to the Fourier-series representation of `ce_m`/`se_m`,
which is equivalent to Parseval's theorem and achieves machine-precision accuracy for
sufficiently large `Nη`.

**Radial integrals (1D numerical):** `I_U` and `I_U′` involve the Modified Mathieu
functions `Ce_m`/`Se_m` over ξ ∈ [0, ξ₀], where ξ₀ = acosh(a/√(a²−b²)). Unlike the
Bessel-based case (where the Lommel integral provides a closed form), no analytic
primitive is known for Modified Mathieu functions, so these are computed by adaptive
Simpson quadrature.

# Arguments
- `a`, `b`: Semi-major and semi-minor axes
- `m`: Azimuthal mode index
- `q`: Mathieu parameter q = (kc·ρ)²/4
- `coeff`: Fourier coefficients of the angular Mathieu function
- `kc`: Cutoff wavenumber
- `β`: Phase constant
- `f`: Frequency in Hz
- `μᵣ`: Relative permeability
- `εᵣ`: Relative permittivity

# Keyword Arguments
- `pol`: `:TE` or `:TM` (default `:TE`)
- `even`: `true` for even (c) modes, `false` for odd (s) modes (default `true`)
- `Nη`: Number of angular quadrature points for `A_v`/`A_vp` (default 4096)
- `rtol`, `atol`, `maxdepth`: Tolerance and depth for the adaptive radial Simpson integrator
"""
function elliptic_modal_power_1d(a, b, m::Int, q, coeff, kc, β, f, μᵣ, εᵣ;
                                            pol=:TE, even::Bool=true,
                                            Nη::Int=4096,
                                            rtol=1e-10, atol=0.0, maxdepth=25)

    ξ0 = elliptic_xi0(a,b)
    _μₒ = 4π*1e-7
    _εₒ = 8.8541878128e-12
    ω = 2π*f
    μ = μᵣ*_μₒ
    ε = εᵣ*_εₒ

    A_v, A_vp = Av_Avp_from_kernels(m, q, coeff; even=even, Nη=Nη)

    radial = even ? Mce_kernel : Mse_kernel

    fU(ξ) = begin
        U, _ = radial(1, m, coeff, q, ξ)
        abs2(U)
    end
    fUp(ξ) = begin
        _, Uξ = radial(1, m, coeff, q, ξ)
        abs2(Uξ)
    end

    Iu  = quad_asimpson(fU,  0.0, ξ0; rtol=rtol, atol=atol, maxdepth=maxdepth)
    Iup = quad_asimpson(fUp, 0.0, ξ0; rtol=rtol, atol=atol, maxdepth=maxdepth)

    I = A_v*Iup + A_vp*Iu

    pref = pol == :TE ? (ω*μ*β)/(kc^4) :
           pol == :TM ? (ω*ε*β)/(kc^4) :
           error("pol must be :TE or :TM")

    return 0.5 * pref * I
end

"""
    te_normalization_ewg(a, b, m, even, kc, β, f, μᵣ, εᵣ; kwargs...)

Normalization factor for TE modes in an elliptical waveguide to achieve unit power.

Returns `√(2/P)`, where `P` is the time-averaged power computed by
[`elliptic_modal_power_1d`]. The Mathieu parameter `q` and
Fourier coefficients are derived internally from `kc`, `m`, and `even`.

The power integral decomposes as:

    P = (ωμβ / 2kc⁴) · (A_v · I_U′ + A_vp · I_U)

- **Angular part** (`A_v`, `A_vp`): quasi-analytical via the spectral trapezoidal rule
  on the Fourier-series of `ce_m`/`se_m` (Parseval-equivalent, see
  [`Av_Avp_from_kernels`]).
- **Radial part** (`I_U`, `I_U′`): 1D adaptive numerical integration of the Modified
  Mathieu functions over ξ ∈ [0, ξ₀] (no closed-form Lommel analog exists).

# Arguments
- `a`, `b`: Semi-major and semi-minor axes
- `m`: Azimuthal mode index
- `even`: `true` for even (c) modes, `false` for odd (s) modes
- `kc`: Cutoff wavenumber (e.g. from `kc_ewg`)
- `β`: Phase constant (e.g. from `phase_constant`)
- `f`: Frequency in Hz
- `μᵣ`: Relative permeability
- `εᵣ`: Relative permittivity

# Keyword Arguments
Forwarded to `elliptic_modal_power_1d` (`Nη`, `rtol`, `atol`, `maxdepth`).
"""
function te_normalization_ewg(a, b, m::Int, even::Bool, kc, β, f, μᵣ, εᵣ; kwargs...)
    ρ = sqrt(a^2 - b^2)
    q = (kc * ρ)^2 / 4
    c = even ? MathieuCharA(m, q) : MathieuCharB(m, q)
    coeff = even ? mathieu_a_coeff(m, q, c, 100) : mathieu_b_coeff(m, q, c, 100)
    P = elliptic_modal_power_1d(a, b, m, q, coeff, kc, β, f, μᵣ, εᵣ;
                                           pol=:TE, even=even, kwargs...)
    return sqrt(2 / P)
end

"""
    tm_normalization_ewg(a, b, m, even, kc, β, f, μᵣ, εᵣ; kwargs...)

Normalization factor for TM modes in an elliptical waveguide to achieve unit power.

Returns `√(2/P)`, where `P` is the time-averaged power computed by
[`elliptic_modal_power_1d`]. The Mathieu parameter `q` and
Fourier coefficients are derived internally from `kc`, `m`, and `even`.

The power integral decomposes as:

    P = (ωεβ / 2kc⁴) · (A_v · I_U′ + A_vp · I_U)

- **Angular part** (`A_v`, `A_vp`): quasi-analytical via the spectral trapezoidal rule
  on the Fourier-series of `ce_m`/`se_m` (Parseval-equivalent, see
  [`Av_Avp_from_kernels`]).
- **Radial part** (`I_U`, `I_U′`): 1D adaptive numerical integration of the Modified
  Mathieu functions over ξ ∈ [0, ξ₀] (no closed-form Lommel analog exists).

# Arguments
- `a`, `b`: Semi-major and semi-minor axes
- `m`: Azimuthal mode index
- `even`: `true` for even (c) modes, `false` for odd (s) modes
- `kc`: Cutoff wavenumber (e.g. from `kc_ewg`)
- `β`: Phase constant (e.g. from `phase_constant`)
- `f`: Frequency in Hz
- `μᵣ`: Relative permeability
- `εᵣ`: Relative permittivity

# Keyword Arguments
Forwarded to `elliptic_modal_power_1d` (`Nη`, `rtol`, `atol`, `maxdepth`).
"""
function tm_normalization_ewg(a, b, m::Int, even::Bool, kc, β, f, μᵣ, εᵣ; kwargs...)
    ρ = sqrt(a^2 - b^2)
    q = (kc * ρ)^2 / 4
    c = even ? MathieuCharA(m, q) : MathieuCharB(m, q)
    coeff = even ? mathieu_a_coeff(m, q, c, 100) : mathieu_b_coeff(m, q, c, 100)
    P = elliptic_modal_power_1d(a, b, m, q, coeff, kc, β, f, μᵣ, εᵣ;
                                           pol=:TM, even=even, kwargs...)
    return sqrt(2 / P)
end

# --- Simpson adaptativo (compacto) ---
@inline function _simpson(f, a, b)
    c = (a + b)/2
    h = b - a
    (h/6) * (f(a) + 4*f(c) + f(b))
end

function _adapt_simpson(f, a, b, fa, fb, fc, S, tol, depth)
    c  = (a + b)/2
    d  = (a + c)/2
    e  = (c + b)/2
    fd = f(d)
    fe = f(e)
    Sleft  = (c - a)/6 * (fa + 4*fd + fc)
    Sright = (b - c)/6 * (fc + 4*fe + fb)
    S2 = Sleft + Sright
    if depth <= 0 || abs(S2 - S) ≤ 15*tol
        return S2 + (S2 - S)/15
    end
    return _adapt_simpson(f, a, c, fa, fc, fd, Sleft,  tol/2, depth-1) +
           _adapt_simpson(f, c, b, fc, fb, fe, Sright, tol/2, depth-1)
end

function quad_asimpson(f, a, b; rtol=1e-10, atol=0.0, maxdepth=20)
    fa = f(a)
    fb = f(b)
    c  = (a + b)/2
    fc = f(c)
    S  = (b - a)/6 * (fa + 4*fc + fb)
    tol = max(atol, rtol*abs(S))
    _adapt_simpson(f, a, b, fa, fb, fc, S, tol, maxdepth)
end