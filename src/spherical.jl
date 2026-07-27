#TODO:     basis = SphericalHarmonics(lmax) is type unstable
# I put everything inside a function barrier and perfomance is good. In any case, it would be good if it were type stable.
#TODO: Define M and N wave vectors, more common on spherical coordinates.

using SpheriCart, StaticArrays

"""
    kc_sph()

Returns 0.0
"""
kc_sph() = 0.0

# spherical bessel functions
@inline jm(m, x) = sphericalbesselj(m, x)
@inline ym(m, x) = sphericalbessely(m, x)
@inline h1m(m, x) = sphericalbesselj(m, x) + im * sphericalbessely(m, x)
@inline h2m(m, x) = sphericalbesselj(m, x) - im * sphericalbessely(m, x)
sph_jm(m, x) = x * jm(m, x)
sph_ym(m, x) = x * ym(m, x)
sph_h1m(m, x) = x * h1m(m, x)
sph_h2m(m, x) = x * h2m(m, x)


"""
    to_svector(x, y, z)

Transform a coordinate or coordinate vectors into a vector of SVector{3,T} 
so that the user does not need to import StaticArrays to calculate problems in 
spherical coordinates with SpheriCart.jl.
"""
function to_svector(x, y, z)
    return SVector(x, y, z)
end

function to_svector(x::AbstractArray{T, N}, y::AbstractArray{T, N}, z::AbstractArray{T, N}) where {T, N}
    coords = similar(x, SVector{3, T})
    for i in eachindex(coords)
        coords[i] = SVector(x[i], y[i], z[i])
    end
    return coords
end

function to_svector(x, y)
    return SVector(x, y)
end

function to_svector(x::AbstractArray{T, N}, y::AbstractArray{T, N}) where {T, N}
    coords = similar(x, SVector{2, T})
    for i in eachindex(coords)
        coords[i] = SVector(x[i], y[i])
    end
    return coords
end

function cart2sph(x, y, z)
    r = hypot(x, y, z)
    ϕ = atan(y, x) 
    θ = atan(hypot(x, y), z)
    return (r, θ, ϕ) 
end

function sph2cart(r, θ, ϕ)
    x = r * sin(θ) * cos(ϕ) 
    y = r * sin(θ) * sin(ϕ) 
    z = r * cos(θ)
    return (x, y, z)
end

function grad_cart2sph(gx, gy, gz, θ, φ)
    gθ     =  cos(θ)*cos(φ)*gx + cos(θ)*sin(φ)*gy - sin(θ)*gz
    gφ     = -sin(φ)*gx        + cos(φ)*gy
    return gθ, gφ
end


"""
    sph_bessel_with_derivatives(f::F, l, x) where F

Return a tuple with f_l(k*r), f'_l(k*r) and (R``+k^2R)  -> l*(l+1)/r^2 * f_l(k*r)

Derivatives of the special spherical bessel functions: `ẑ(r) = k*r * zm(m, k*r)`.

"""
function sph_bessel_with_derivatives(f::F, l, r, k) where F <: Function
    x = r * k
    h  = f(l, x)
    h′ = (l+1)/r * h - k * f(l+1, x)  
    h″_k = l*(l+1)/r^2 * h
    return (h, h′, h″_k)
end



sph_jm_with_derivatives(l, r, k) = sph_bessel_with_derivatives(sph_jm, l, r, k) 
sph_ym_with_derivatives(l, r, k) = sph_bessel_with_derivatives(sph_ym, l, r, k) 
sph_h1m_with_derivatives(l, r, k) = sph_bessel_with_derivatives(sph_h1m, l, r, k) 
sph_h2m_with_derivatives(l, r, k) = sph_bessel_with_derivatives(sph_h2m, l, r, k) 


function sph_modal_f(θ, φ, rs, ylm, ylm_p, k)

    R, R′, R″_k = rs

    gx, gy, gz = ylm_p
    gθ, gφ = grad_cart2sph(gx, gy, gz, θ, φ)

    ψ   = R * ylm
    ∂ψθ = R * gθ              
    ∂ψφ = R * gφ  

    ∂²ψᵣθ = R′ * gθ     # ∂²Fᵣ/∂r∂θ
    ∂²ψᵣφ = R′ * gφ     # ∂²Fᵣ/∂r∂φ
    ∂²ψᵣᵣ = R″_k * ylm  # ∂²Fᵣ/∂r² + k²Fᵣ

    return (ψ, ∂ψθ, ∂ψφ, ∂²ψᵣθ, ∂²ψᵣφ ,∂²ψᵣᵣ)
end

"""
    te_sph_fields(r, θ, ϕ, rs, ylm, ylm_p, k, μᵣ, εᵣ)

Compute TE spherical mode fields at a single point `(r, θ, ϕ)` given pre-computed radial
data `rs = (R, R′, R″_k)`, spherical harmonic value `ylm`, its gradient `ylm_p`, and
medium parameters. Returns `(Eᵣ, Eθ, Eϕ, Hᵣ, Hθ, Hϕ)`.
"""
function te_sph_fields(r, θ, ϕ, rs, ylm, ylm_p, k, μᵣ, εᵣ)

    μ = μᵣ * _μₒ
    ε = εᵣ * _εₒ
    ω = k / sqrt(μ * ε)

    ψ, ∂ψθ, ∂ψφ, ∂²ψᵣθ, ∂²ψᵣφ ,∂²ψᵣᵣ = sph_modal_f(θ, ϕ, rs, ylm, ylm_p, k)

    E_θ = -1/ε * ∂ψφ
    E_ϕ =  1/ε * ∂ψθ
    E_r = zero(E_θ)
    
    H_r = 1/(im*ω*μ*ε) * ∂²ψᵣᵣ   # k^2 * ψ is already on ∂²ψᵣᵣ  1/(im*ω*μ*ε) * (∂²ψᵣᵣ + k^2 * ψ)
    H_θ = 1/(im*ω*μ*ε) * ∂²ψᵣθ
    H_ϕ = 1/(im*ω*μ*ε) * ∂²ψᵣφ

    return (E_r, E_θ, E_ϕ, H_r, H_θ, H_ϕ)
end


function te_sph_fields(r_vec::SVector{3, T}, rs, ylm, ylm_p, k, μᵣ, εᵣ) where T
    x, y, z = r_vec
    r, θ, ϕ = cart2sph(x, y, z)
    return te_sph_fields(r, θ, ϕ, rs, ylm, ylm_p, k, μᵣ, εᵣ)
end


function te_sph_fields(r_vec::Vector{SVector{3, T}}, rs, ylm, ylm_p, k, μᵣ, εᵣ) where T
    fields = similar(r_vec, NTuple{6, Complex{T}})
    for idx in eachindex(r_vec)
        x, y, z = r_vec[idx]
        r, θ, ϕ = cart2sph(x, y, z)
        fields[idx] = te_sph_fields(r, θ, ϕ, rs[idx], ylm[idx], ylm_p[idx], k, μᵣ, εᵣ)
    end
    return fields
end


"""
    te_sph_fields_lmax(r_vec, lmax, f, μᵣ, εᵣ, radial=4)

Compute TE spherical mode fields for all `(l, m)` up to `lmax`.

Returns a matrix `A[p, i]` of `NTuple{6, ComplexF64}` with `(Eᵣ, Eθ, Eϕ, Hᵣ, Hθ, Hϕ)` in
spherical components. The mode index follows `i = l² + l + m + 1`.

# Arguments
- `r_vec`: vector of `SVector{3, T}` with point coordinates in Cartesian form.
- `lmax`: maximum degree `l`. Total number of modes is `(lmax+1)²`.
- `f`: frequency in Hz.
- `μᵣ`: relative permeability.
- `εᵣ`: relative permittivity.
- `radial`: radial function type (1=j, 2=y, 3=h₁ incoming, 4=h₂ outgoing). Default: `4`.
"""
function te_sph_fields_lmax(r_vec, lmax::Int, f, μᵣ, εᵣ, radial::Int = 4)

    basis = SphericalHarmonics(lmax, normalisation = :sphericart)
    A = Array{NTuple{6, ComplexF64}, 2}(undef, length(r_vec), (lmax+1)^2)
    Rs = Matrix{NTuple{3, ComplexF64}}(undef, length(r_vec), lmax+1)
    te_sph_fields_lmax!(A, Rs, r_vec, basis, lmax, f, μᵣ, εᵣ, radial)

    return A
end


function te_sph_fields_lmax!(A, Rs, r_vec, basis, lmax::Int, f, μᵣ, εᵣ, radial::Int)

    sph_coords = map(x->cart2sph(x[1], x[2], x[3]),r_vec)

    k = wavenumber(f, μᵣ, εᵣ)
    for l in 0:lmax
        for i in eachindex(r_vec)
            x, y, z = r_vec[i]
            r = hypot(x, y, z)
            if radial == 3
                R, R´, R´´ = sph_h1m_with_derivatives(l, r, k)
            else
                R, R´, R´´ = sph_h2m_with_derivatives(l, r, k)
            end
            Rs[i, l+1] = (R, R´, R´´)
        end
    end

    Ylm, ∇Ylm = compute_with_gradients(basis, r_vec)

    for l in 0:lmax
        for m in -l:l
            #idx = SpheriCart.lm2idx(l, m) # not public
            idx = l^2 + l + m + 1
            for ijk in eachindex(r_vec)
                r, θ, ϕ = sph_coords[ijk]
                A[ijk, idx] = te_sph_fields(r, θ, ϕ, Rs[ijk, l+1], Ylm[ijk, idx], ∇Ylm[ijk, idx], k, μᵣ, εᵣ)
            end
        end
    end
    return nothing
end


"""
    tm_sph_fields(r, θ, ϕ, rs, ylm, ylm_p, k, μᵣ, εᵣ)

Compute TM spherical mode fields at a single point `(r, θ, ϕ)` given pre-computed radial
data `rs = (R, R′, R″_k)`, spherical harmonic value `ylm`, its gradient `ylm_p`, and
medium parameters. Returns `(Eᵣ, Eθ, Eϕ, Hᵣ, Hθ, Hϕ)`.
"""
function tm_sph_fields(r, θ, ϕ, rs, ylm, ylm_p, k, μᵣ, εᵣ)
    μ = μᵣ * _μₒ
    ε = εᵣ * _εₒ
    ω = k / sqrt(μ * ε)
    ψ, ∂ψθ, ∂ψφ, ∂²ψᵣθ, ∂²ψᵣφ, ∂²ψᵣᵣ = sph_modal_f(θ, ϕ, rs, ylm, ylm_p, k)
    
    E_r = 1/(im*ω*μ*ε) * ∂²ψᵣᵣ # k^2 * ψ is already on ∂²ψᵣᵣ /// 1/(im*ω*μ*ε) * (∂²ψᵣᵣ + k^2 * ψ)
    E_θ = 1/(im*ω*μ*ε) * ∂²ψᵣθ
    E_ϕ = 1/(im*ω*μ*ε) * ∂²ψᵣφ
    
    H_r = zero(E_θ)  
    H_θ =  1/μ * ∂ψφ
    H_ϕ = -1/μ * ∂ψθ

    return (E_r, E_θ, E_ϕ, H_r, H_θ, H_ϕ)
end


"""
    tm_sph_fields_lmax(r_vec, lmax, f, μᵣ, εᵣ, radial=4)

Compute TM spherical mode fields for all `(l, m)` up to `lmax`.

Returns a matrix `A[p, i]` of `NTuple{6, ComplexF64}` with `(Eᵣ, Eθ, Eϕ, Hᵣ, Hθ, Hϕ)` in
spherical components. The mode index follows `i = l² + l + m + 1`.

# Arguments
- `r_vec`: vector of `SVector{3, T}` with point coordinates in Cartesian form.
- `lmax`: maximum degree `l`. Total number of modes is `(lmax+1)²`.
- `f`: frequency in Hz.
- `μᵣ`: relative permeability.
- `εᵣ`: relative permittivity.
- `radial`: radial function type (1=j, 2=y, 3=h₁ incoming, 4=h₂ outgoing). Default: `4`.
"""
function tm_sph_fields_lmax(r_vec, lmax::Int, f, μᵣ, εᵣ, radial::Int = 4)

    basis = SphericalHarmonics(lmax, normalisation = :sphericart)
    A = Array{NTuple{6, ComplexF64}, 2}(undef, length(r_vec), (lmax+1)^2)
    Rs = Matrix{NTuple{3, ComplexF64}}(undef, length(r_vec), lmax+1)
    tm_sph_fields_lmax!(A, Rs, r_vec, basis, lmax, f, μᵣ, εᵣ, radial)

    return A
end


function tm_sph_fields_lmax!(A, Rs, r_vec, basis, lmax::Int, f, μᵣ, εᵣ, radial::Int)

    sph_coords = map(x->cart2sph(x[1], x[2], x[3]),r_vec)

    k = wavenumber(f, μᵣ, εᵣ)
    for l in 0:lmax
        for i in eachindex(r_vec)
            x, y, z = r_vec[i]
            r = hypot(x, y, z)
            if radial == 3
                R, R´, R´´ = sph_h1m_with_derivatives(l, r, k)
            else
                R, R´, R´´ = sph_h2m_with_derivatives(l, r, k)
            end
            Rs[i, l+1] = (R, R´, R´´)
        end
    end

    Ylm, ∇Ylm = compute_with_gradients(basis, r_vec)

    for l in 0:lmax
        for m in -l:l
            #idx = SpheriCart.lm2idx(l, m) # not public
            idx = l^2 + l + m + 1
            for ijk in eachindex(r_vec)
                r, θ, ϕ = sph_coords[ijk]
                A[ijk, idx] = tm_sph_fields(r, θ, ϕ, Rs[ijk, l+1], Ylm[ijk, idx], ∇Ylm[ijk, idx], k, μᵣ, εᵣ)
            end
        end
    end
    return nothing
end

function tm_sph_fields(r_vec, rs, ylm, ylm_p, k, μᵣ, εᵣ)
    x, y, z = r_vec
    r, θ, ϕ = cart2sph(x, y, z)
    return tm_sph_fields(r, θ, ϕ, rs, ylm, ylm_p, k, μᵣ, εᵣ)
end


function tm_sph_fields(r_vec::Vector{SVector{3, T}}, rs, ylm, ylm_p, k, μᵣ, εᵣ) where T
    fields = similar(r_vec, NTuple{6, Complex{T}})
    for idx in eachindex(r_vec)
        x, y, z = r_vec[idx]
        r, θ, ϕ = cart2sph(x, y, z)
        fields[idx] = tm_sph_fields(r, θ, ϕ, rs, ylm[idx], ylm_p[idx], k, μᵣ, εᵣ)
    end
    return fields
end

"""
    find_lmax_n(N)

Return the minimum `lmax` that produces N modes (using TE and TM modes).
"""
find_lmax_n(N) = ceil(Int, sqrt(N/2) - 1)

"""
    first_n_modes_sph(N)

Returns the first `N` spherical modes.
"""
function first_n_modes_sph(N)
    modes = []
    lmax = find_lmax_n(N)
    push!(modes, (:TE, 0, 0))
    push!(modes, (:TM, 0, 0))
    
    for l in 1:lmax
        ms_sorted = sort(-l:l, by=abs)  # [0, -1, 1, -2, 2, ...]
        
        for m in ms_sorted
            push!(modes, (:TE, l, m))
            push!(modes, (:TM, l, m))
        end
    end
    
    return modes[1:N]
end


function metric_and_unit_spherical(r, θ, ϕ)
    st, ct = sincos(θ)
    sp, cp = sincos(ϕ)
    
    h_r = 1.0
    h_θ = r
    h_ϕ = r * st
    
    e_r_x = st * cp
    e_r_y = st * sp  
    e_r_z = ct
    
    e_θ_x = ct * cp
    e_θ_y = ct * sp
    e_θ_z = -st
    
    e_ϕ_x = -sp
    e_ϕ_y = cp
    e_ϕ_z = 0.0
    
    return h_r, h_θ, h_ϕ, e_r_x, e_r_y, e_r_z, e_θ_x, e_θ_y, e_θ_z, e_ϕ_x, e_ϕ_y, e_ϕ_z
end


"""
    spherical_to_cartesian_fields(E_r, E_θ, E_ϕ, H_r, H_θ, H_ϕ, θ, ϕ)

Convert electromagnetic field components from spherical `(r, θ, ϕ)` to Cartesian `(x, y, z)`.
Returns `(Eₓ, Eᵧ, E_z, Hₓ, Hᵧ, H_z)`.
"""
function spherical_to_cartesian_fields(E_r, E_θ, E_ϕ, H_r, H_θ, H_ϕ, θ, ϕ)

    st, ct = sincos(θ)
    sp, cp = sincos(ϕ)
    
    # Transformación directa
    E_x = E_r * st*cp + E_θ * ct*cp - E_ϕ * sp
    E_y = E_r * st*sp + E_θ * ct*sp + E_ϕ * cp
    E_z = E_r * ct     - E_θ * st
    
    H_x = H_r * st*cp + H_θ * ct*cp - H_ϕ * sp  
    H_y = H_r * st*sp + H_θ * ct*sp + H_ϕ * cp
    H_z = H_r * ct     - H_θ * st
    
    return (E_x, E_y, E_z, H_x, H_y, H_z)
end


"""
    tm_normalization_sph(l, r, k, radial, μᵣ, εᵣ)

Normalization factor for TM spherical modes to achieve unit power.

Derived by integrating the Poynting vector over a sphere of radius `r`.

# Arguments
- `l`: degree of the mode.
- `r`: evaluation radius.
- `k`: wavenumber.
- `radial`: radial function type (3=h₁, 4=h₂).
- `μᵣ`: relative permeability.
- `εᵣ`: relative permittivity.
"""
function tm_normalization_sph(l, r, k, radial::Int, μᵣ, εᵣ)
    ε = εᵣ * _εₒ
    μ = μᵣ * _μₒ
    c = 1 / sqrt(μ * ε)
    ω = k * c
    R, Rp, _ = radial == 3 ? sph_h1m_with_derivatives(l, r, k) : sph_h2m_with_derivatives(l, r, k)
    P_unnormalized = abs( imag(R * conj(Rp)) / (2 * ω * μ^2 * ε) * l * (l+1))
    norm_factor = sqrt(1 / P_unnormalized)
    return norm_factor
end

    
"""
    te_normalization_sph(l, r, k, radial, μᵣ, εᵣ)

Normalization factor for TE spherical modes to achieve unit power.

Derived by integrating the Poynting vector over a sphere of radius `r`.

# Arguments
- `l`: degree of the mode.
- `r`: evaluation radius.
- `k`: wavenumber.
- `radial`: radial function type (3=h₁, 4=h₂).
- `μᵣ`: relative permeability.
- `εᵣ`: relative permittivity.
"""
function te_normalization_sph(l, r, k, radial::Int, μᵣ, εᵣ)
    ε = εᵣ * _εₒ
    μ = μᵣ * _μₒ
    c = 1 / sqrt(μ * ε)
    ω = k * c
    R, Rp, _ = radial == 3 ? sph_h1m_with_derivatives(l, r, k) : sph_h2m_with_derivatives(l, r, k)
    P_unnormalized = abs( imag(R * conj(Rp)) / (2 * ω * μ * ε^2) * l * (l+1))
    norm_factor = sqrt(1 / P_unnormalized)
    return norm_factor
end


"""
    m_normalization_sph(l, r, k, radial, μᵣ, εᵣ)

**Basis normalization**, not power normalization.

Normalization factor for the spherical vector wave function **Mₗₘ**
based on surface L² normalization on a sphere:

    ∫ |Mₗₘ|² dA = 1

after scaling by this factor.

This normalization only fixes the amplitude of **M** as a basis function.
It does **not** correspond to unit electromagnetic power.

In the spherical vector-wave-function formalism, physical power is associated
with the bilinear flux pairing between **M** and **N**:

    P ∝ ∫ Im{ Mₜ × Nₜ* } · r̂ dA

The complementary scaling between **M** and **N** is therefore set through
`n_normalization_sph`, which fixes the M–N flux normalization.

# Arguments
- `l`: degree of the mode.
- `r`: evaluation radius.
- `k`: wavenumber.
- `radial`: radial function type (3=h₁, 4=h₂).
- `μᵣ`: relative permeability.
- `εᵣ`: relative permittivity.
"""
function m_normalization_sph(l, r, k, radial::Int, μᵣ, εᵣ)
    R, _, _ = radial == 3 ? sph_h1m_with_derivatives(l, r, k) : sph_h2m_with_derivatives(l, r, k)
    norm_factor = sqrt(1 / (l*(l+1) * abs2(R)))
    return norm_factor
end

"""
    n_normalization_sph(l, r, k, radial, μᵣ, εᵣ)

**Basis normalization**, not power normalization.

Normalization factor for the spherical vector wave function **Nₗₘ**
defined so that, together with the scaling returned by
`m_normalization_sph`, the bilinear M–N flux pairing on a sphere is
unit-normalized.

This function does **not** normalize **N** independently in an `L²` sense.
Instead, it fixes the relative scaling between **M** and **N** through the
radial flux invariant:

    ∫ 1/2 * Im{ Mₜ × Nₜ* } · r̂ dA = 1

after applying both normalization factors.

In practice:

- `m_normalization_sph` first fixes the amplitude of **M** through a
  surface `L²` normalization (`∫ |M|² dA = 1`);
- `n_normalization_sph` then determines the complementary scaling of **N**
  so that the bilinear M–N pairing becomes unitary.

Therefore, this normalization **depends on the convention used for**
`m_normalization_sph`. The two functions must be used together.

This is a **basis normalization convention** for the spherical
vector-wave-function formalism. It is **not** equivalent to unit
electromagnetic power normalization of a physical TE/TM mode.

For physical unit-power scaling, use the TE/TM normalization functions instead.

# Arguments
- `l`: degree of the mode.
- `r`: evaluation radius.
- `k`: wavenumber.
- `radial`: radial function type (3=h₁, 4=h₂).
- `μᵣ`: relative permeability.
- `εᵣ`: relative permittivity.
"""
function n_normalization_sph(l, r, k, radial::Int, μᵣ, εᵣ)
    R, Rp, _ = radial == 3 ? sph_h1m_with_derivatives(l, r, k) : sph_h2m_with_derivatives(l, r, k)
    Nm = m_normalization_sph(l, r, k, radial, μᵣ, εᵣ)
    Pmn_unnormalized = abs(imag(R * conj(Rp)) / (2 * k) * l * (l+1))
    norm_factor = 1 / (Pmn_unnormalized * Nm)
    return norm_factor
end



"""
    mn_vectors_sph(r, θ, ϕ, rs, ylm, ylm_p, k, μᵣ, εᵣ)

Computes the **M_{lm}** and **N_{lm}** spherical wave vectors.
"""
function mn_vectors_sph(r, θ, ϕ, rs, ylm, ylm_p, k, μᵣ, εᵣ)

    ψ, ∂ψθ, ∂ψφ, ∂²ψᵣθ, ∂²ψᵣφ ,∂²ψᵣᵣ = sph_modal_f(θ, ϕ, rs, ylm, ylm_p, k)
 
    mφ = ∂ψφ
    mθ = -∂ψθ
    mr = zero(mφ)
    nr = 1/k * ∂²ψᵣᵣ
    nθ = 1/k * ∂²ψᵣφ
    nφ = 1/k * ∂²ψᵣθ

    return (mr, mθ, mφ, nr, nθ, nφ)
end

"""
    te_from_mn_sph(r, θ, ϕ, rs, ylm, ylm_p, l, k, radial, μᵣ, εᵣ; normalize=true)

Build TE spherical fields from the `M/N` vector basis and return
`(Eᵣ, Eθ, Eϕ, Hᵣ, Hθ, Hϕ)`.

When `normalize=true` (default), fields are scaled with
`te_normalization_sph` to unit power.
"""
function te_from_mn_sph(r, θ, ϕ, rs, ylm, ylm_p, l, k, radial::Int, μᵣ, εᵣ; normalize::Bool = true)
    μ = μᵣ * _μₒ
    ε = εᵣ * _εₒ
    ω = k / sqrt(μ * ε)

    Mᵣ, Mθ, Mϕ, Nᵣ, Nθ, Nϕ = mn_vectors_sph(r, θ, ϕ, rs, ylm, ylm_p, k, μᵣ, εᵣ)

    # Match the same component convention already used by `te_sph_fields`.
    Eᵣ = zero(Mᵣ)
    Eθ = -1/ε * Mϕ
    Eϕ = -1/ε * Mθ

    sN = k / (im * ω * μ * ε)
    Hᵣ = sN * Nᵣ
    Hθ = sN * Nϕ
    Hϕ = sN * Nθ

    A = normalize ? te_normalization_sph(l, r, k, radial, μᵣ, εᵣ) : 1.0
    return (A * Eᵣ, A * Eθ, A * Eϕ, A * Hᵣ, A * Hθ, A * Hϕ)
end

"""
    tm_from_mn_sph(r, θ, ϕ, rs, ylm, ylm_p, l, k, radial, μᵣ, εᵣ; normalize=true)

Build TM spherical fields from the `M/N` vector basis and return
`(Eᵣ, Eθ, Eϕ, Hᵣ, Hθ, Hϕ)`.

When `normalize=true` (default), fields are scaled with
`tm_normalization_sph` to unit power.
"""
function tm_from_mn_sph(r, θ, ϕ, rs, ylm, ylm_p, l, k, radial::Int, μᵣ, εᵣ; normalize::Bool = true)
    μ = μᵣ * _μₒ
    ε = εᵣ * _εₒ
    ω = k / sqrt(μ * ε)

    Mᵣ, Mθ, Mϕ, Nᵣ, Nθ, Nϕ = mn_vectors_sph(r, θ, ϕ, rs, ylm, ylm_p, k, μᵣ, εᵣ)

    sN = k / (im * ω * μ * ε)
    Eᵣ = sN * Nᵣ
    Eθ = sN * Nϕ
    Eϕ = sN * Nθ

    Hᵣ = zero(Mᵣ)
    Hθ = 1/μ * Mϕ
    Hϕ = 1/μ * Mθ

    A = normalize ? tm_normalization_sph(l, r, k, radial, μᵣ, εᵣ) : 1.0
    return (A * Eᵣ, A * Eθ, A * Eϕ, A * Hᵣ, A * Hθ, A * Hϕ)
end


"""
    mn_sph_vectors_lmax(r_vec, lmax, f, μᵣ, εᵣ, radial=4)

Compute normalized **M** and **N** spherical wave vectors for all `(l, m)` up to `lmax`.

Returns a matrix `A[p, i]` of `NTuple{6, ComplexF64}` with `(Mᵣ, Mθ, Mϕ, Nᵣ, Nθ, Nϕ)` in
spherical components. Modes are indexed by `i = l² + l + m + 1`.

# Arguments
- `r_vec`: vector of `SVector{3, T}` with point coordinates in Cartesian form.
- `lmax`: maximum degree `l`. Total number of modes is `(lmax+1)²`.
- `f`: frequency in Hz.
- `μᵣ`: relative permeability.
- `εᵣ`: relative permittivity.
- `radial`: radial function type (1=j, 2=y, 3=h₁ incoming, 4=h₂ outgoing). Default: `4`.
"""
function mn_sph_vectors_lmax(r_vec, lmax::Int, f, μᵣ, εᵣ, radial::Int = 4)

    basis = SphericalHarmonics(lmax, normalisation = :sphericart)
    A = Array{NTuple{6, ComplexF64}, 2}(undef, length(r_vec), (lmax+1)^2)
    Rs = Matrix{NTuple{3, ComplexF64}}(undef, length(r_vec), lmax+1)
    mn_sph_vectors_lmax!(A, Rs, r_vec, basis, lmax, f, μᵣ, εᵣ, radial)

    return A
end

"""
    m_sph_vectors_lmax(r_vec, lmax, f, μᵣ, εᵣ, radial=4)

Compute normalized **M** spherical wave vectors for all `(l, m)` up to `lmax`.
Returns a matrix of `NTuple{3, ComplexF64}` with `(Mᵣ, Mθ, Mϕ)`. See `mn_sph_vectors_lmax`.
"""
function m_sph_vectors_lmax(r_vec, lmax::Int, f, μᵣ, εᵣ, radial::Int = 4)
    B = mn_sph_vectors_lmax(r_vec, lmax, f, μᵣ, εᵣ, radial)
    return [(b[1], b[2], b[3]) for b in B]
end

"""
    n_sph_vectors_lmax(r_vec, lmax, f, μᵣ, εᵣ, radial=4)

Compute normalized **N** spherical wave vectors for all `(l, m)` up to `lmax`.
Returns a matrix of `NTuple{3, ComplexF64}` with `(Nᵣ, Nθ, Nϕ)`. See `mn_sph_vectors_lmax`.
"""
function n_sph_vectors_lmax(r_vec, lmax::Int, f, μᵣ, εᵣ, radial::Int = 4)
    B = mn_sph_vectors_lmax(r_vec, lmax, f, μᵣ, εᵣ, radial)
    return [(b[4], b[5], b[6]) for b in B]
end

function mn_sph_vectors_lmax!(A, Rs, r_vec, basis, lmax::Int, f, μᵣ, εᵣ, radial::Int)

    sph_coords = map(x->cart2sph(x[1], x[2], x[3]),r_vec)
    xi, yi, zi = r_vec[1]
    ri = hypot(xi, yi, zi)
    k = wavenumber(f, μᵣ, εᵣ)
    for l in 0:lmax
        for i in eachindex(r_vec)
            x, y, z = r_vec[i]
            r = hypot(x, y, z)
            if radial == 3
                R, R´, R´´ = sph_h1m_with_derivatives(l, r, k)
            else
                R, R´, R´´ = sph_h2m_with_derivatives(l, r, k)
            end
            Rs[i, l+1] = (R, R´, R´´)
        end
    end

    Ylm, ∇Ylm = compute_with_gradients(basis, r_vec)

    for l in 1:lmax
        factor = m_normalization_sph(l, ri, k, radial, μᵣ, εᵣ)
        for m in -l:l
            #idx = SpheriCart.lm2idx(l, m) # not public
            idx = l^2 + l + m + 1
            for ijk in eachindex(r_vec)
                r, θ, ϕ = sph_coords[ijk]
                A[ijk, idx] = mn_vectors_sph(r, θ, ϕ, Rs[ijk, l+1], Ylm[ijk, idx], ∇Ylm[ijk, idx], k, μᵣ, εᵣ) .* factor
            end
        end
    end
    return nothing
end
