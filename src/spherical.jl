#TODO: Use normalisation = :sphericart -  the same as `:L2`, gives L2-orthonormality, i.e. ∫ |Ylm|² = 1?
#TODO:     basis = SphericalHarmonics(lmax) is type unstable
# I put everything inside a function barrier and perfomance is good. In any case, it would be good if it were type stable.

using SpheriCart, StaticArrays

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
    to_svector3(x, y, z)

Transform a coordinate or coordinate vectors into a vector of SVector{3,T} 
so that the user does not need to import StaticArrays to calculate problems in 
spherical coordinates with SpheriCart.jl.
"""
function to_svector3(x, y, z)
    return SVector(x, y, z)
end

function to_svector3(x::AbstractArray{T, N}, y::AbstractArray{T, N}, z::AbstractArray{T, N}) where {T, N}
    coords = similar(x, SVector{3, T})
    for i in eachindex(coords)
        coords[i] = SVector(x[i], y[i], z[i])
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

Return a tuple with f_l(x), f'_l(x) and (R``+k^2R)  -> l*(l+1)/x^2 * f_l(x)
"""
function sph_bessel_with_derivatives(f::F, l, x) where F <: Function
    h  = f(l, x)
    h′ = (l/x) * h - f(l+1, x)  
    #h″ = ((l*(l+1))/x^2 - 1) * h - (2/x) * h′
    # return (R``+k^2)Fr
    h″ = l*(l+1)/x^2 * h
    return (h, h′, h″)
end

sph_jm_with_derivatives(l, x) = sph_bessel_with_derivatives(sph_jm, l, x) 
sph_ym_with_derivatives(l, x) = sph_bessel_with_derivatives(sph_ym, l, x) 
sph_h1m_with_derivatives(l, x) = sph_bessel_with_derivatives(sph_h1m, l, x) 
sph_h2m_with_derivatives(l, x) = sph_bessel_with_derivatives(sph_h2m, l, x) 


function sph_modal_f(θ, φ, rs, ylm, ylm_p, k)

    R, R′, R″ = rs

    gx, gy, gz = ylm_p
    gθ, gφ = grad_cart2sph(gx, gy, gz, θ, φ)

    ψ   = R * ylm
    ∂ψθ = R * gθ              
    ∂ψφ = R * gφ  

    ∂²ψᵣθ = k * R′ * gθ     # ∂²Fᵣ/∂r∂θ
    ∂²ψᵣφ = k * R′ * gφ     # ∂²Fᵣ/∂r∂φ
    ∂²ψᵣᵣ = R″ * ylm  # ∂²Fᵣ/∂r² + k²Fᵣ

    return (ψ, ∂ψθ, ∂ψφ, ∂²ψᵣθ, ∂²ψᵣφ ,∂²ψᵣᵣ)
end

#TODO: move ω outside
function te_sph_fields(r, θ, ϕ, rs, ylm, ylm_p, k, μᵣ, εᵣ)

    μ = μᵣ * _μₒ
    ε = εᵣ * _εₒ
    ω = k / sqrt(μ * ε)

    ψ, ∂ψθ, ∂ψφ, ∂²ψᵣθ, ∂²ψᵣφ ,∂²ψᵣᵣ = sph_modal_f(θ, ϕ, rs, ylm, ylm_p, k)

    E_θ = -1/ε * ∂ψφ
    E_ϕ =  1/ε * ∂ψθ
    E_r = zero(E_θ)
    
    H_r = 1/(im*ω*μ*ε) * (∂²ψᵣᵣ + k^2 * ψ)
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
    te_sph_fields_lmax(r_vec, lmax, f, μᵣ, εᵣ, incident)

# Arguments
- `r_vec`: vector of `SVector{3, T}` with the position of each point in cartesian coordinates.
- `lmax`: max value of L of the modes TE_{lm}. The total number of modes is (lmax+1)^2.
- `f`: Frequency in Hz
- `μᵣ`: Relative permeability
- `εᵣ`: Relative permittivity
- `incident`: ongoing or outgoing wave (`true` or `false`)

"""
function te_sph_fields_lmax(r_vec, lmax::Int, f, μᵣ, εᵣ, incident)

    basis = SphericalHarmonics(lmax, normalisation = :sphericart)
    A = Array{NTuple{6, ComplexF64}, 2}(undef, length(r_vec), (lmax+1)^2)
    Rs = Matrix{NTuple{3, ComplexF64}}(undef, length(r_vec), lmax+1)
    _te_sph_fields_lmax!(A, Rs, r_vec, basis, lmax, f, μᵣ, εᵣ, incident)
    
    return A
end


function _te_sph_fields_lmax!(A, Rs, r_vec, basis, lmax::Int, f, μᵣ, εᵣ, incident)


    sph_coords = map(x->cart2sph(x[1], x[2], x[3]),r_vec)

    k = wavenumber(f, μᵣ, εᵣ)
    for l in 0:lmax
        for i in eachindex(r_vec) 
            x, y, z = r_vec[i]
            r = hypot(x, y, z)
            R, R´, R´´ =  incident == true ? sph_h1m_with_derivatives(l, k*r) : sph_h2m_with_derivatives(l, k*r)
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

"""
function tm_sph_fields(r, θ, ϕ, rs, ylm, ylm_p, k, μᵣ, εᵣ)
    μ = μᵣ * _μₒ
    ε = εᵣ * _εₒ
    ω = k / sqrt(μ * ε)
    ψ, ∂ψθ, ∂ψφ, ∂²ψᵣθ, ∂²ψᵣφ, ∂²ψᵣᵣ = sph_modal_f(θ, ϕ, rs, ylm, ylm_p, k)
    
    E_r = 1/(im*ω*μ*ε) * (∂²ψᵣᵣ + k^2 * ψ)
    E_θ = 1/(im*ω*μ*ε) * ∂²ψᵣθ
    E_ϕ = 1/(im*ω*μ*ε) * ∂²ψᵣφ
    
    H_r = zero(E_θ)  
    H_θ =  1/μ * ∂ψφ
    H_ϕ = -1/μ * ∂ψθ

    return (E_r, E_θ, E_ϕ, H_r, H_θ, H_ϕ)
end


"""
    tm_sph_fields_lmax(r_vec, lmax, f, μᵣ, εᵣ, incident)

# Arguments
- `r_vec`: vector of `SVector{3, T}` with the position of each point in cartesian coordinates.
- `lmax`: max value of L of the modes TM_{lm}. The total number of modes is (lmax+1)^2.
- `f`: Frequency in Hz
- `μᵣ`: Relative permeability
- `εᵣ`: Relative permittivity
- `incident::Bool`: ongoing or outgoing wave.

"""
function tm_sph_fields_lmax(r_vec, lmax::Int, f, μᵣ, εᵣ, incident)

    basis = SphericalHarmonics(lmax, normalisation = :sphericart)
    A = Array{NTuple{6, ComplexF64}, 2}(undef, length(r_vec), (lmax+1)^2)
    Rs = Matrix{NTuple{3, ComplexF64}}(undef, length(r_vec), lmax+1)
    _tm_sph_fields_lmax!(A, Rs, r_vec, basis, lmax, f, μᵣ, εᵣ, incident)
    
    return A
end


function _tm_sph_fields_lmax!(A, Rs, r_vec, basis, lmax::Int, f, μᵣ, εᵣ, incident)


    sph_coords = map(x->cart2sph(x[1], x[2], x[3]),r_vec)

    k = wavenumber(f, μᵣ, εᵣ)
    for l in 0:lmax
        for i in eachindex(r_vec) 
            x, y, z = r_vec[i]
            r = hypot(x, y, z)
            R, R´, R´´ =  incident == true ? sph_h1m_with_derivatives(l, k*r) : sph_h2m_with_derivatives(l, k*r)
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
    tm_normalization_sph(r, k, incident, μᵣ, εᵣ)

Normalization factor for TM modes to achieve unit power.

The expression can be derived by integrating the Poynting vector over the cross-section of the guide.

# Arguments
- `l`: mode index
- `r`: Radial coordinate
- `k`: propagation constant
- `incident::Bool`: ongoing or outgoing wave 
- `μᵣ`: Relative permeability
- `εᵣ`: Relative permittivity
"""
function tm_normalization_sph(l, r, k, incident, μᵣ, εᵣ)
    ε = εᵣ * _εₒ
    μ = μᵣ * _μₒ
    c = 1 / sqrt(μ * ε)
    ω = k * c
    R, Rp, _ = incident == true ? sph_h1m_with_derivatives(l, k*r) : sph_h2m_with_derivatives(l, k*r)
    P_unnormalized = abs( (k * r^2)/(2 * ω * μ^2 * ε) * real(im * R * Rp) * l * (l+1))
    norm_factor = sqrt(1 / P_unnormalized)
    return norm_factor
end

    
"""
    te_normalization_sph(l, r, k, incident, μᵣ, εᵣ)

Normalization factor for TE modes to achieve unit power.

The expression can be derived by integrating the Poynting vector over the cross-section of the guide.

# Arguments
- `l`: mode index
- `r`: Radial coordinate
- `k`: propagation constant
- `incident::Bool`: ongoing or outgoing wave 
- `μᵣ`: Relative permeability
- `εᵣ`: Relative permittivity
"""
function te_normalization_sph(l, r, k, incident, μᵣ, εᵣ)
    ε = εᵣ * _εₒ
    μ = μᵣ * _μₒ
    c = 1 / sqrt(μ * ε)
    ω = k * c
    R, Rp, _ = incident == true ? sph_h1m_with_derivatives(l, k*r) : sph_h2m_with_derivatives(l, k*r)
    P_unnormalized = abs( (k * r^2)/(2 * ω * μ * ε^2) * real(im * R * Rp) * l * (l+1))
    norm_factor = sqrt(1 / P_unnormalized)
    return norm_factor
end

